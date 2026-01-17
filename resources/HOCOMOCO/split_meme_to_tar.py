#!/usr/bin/env python3
import argparse
import os
import re
import tarfile
import tempfile
from io import BytesIO
from typing import List, Tuple


def read_file(path: str) -> List[str]:
    with open(path, "r") as f:
        return f.readlines()


def split_header_and_motifs(lines: List[str]) -> Tuple[List[str], List[List[str]]]:
    """
    Split a MEME v4 file into:
    - header_lines: everything up to (but not including) the first 'MOTIF' line
    - motifs: list of motif blocks (each starting with 'MOTIF ...')
    """
    header_lines = []
    motifs = []

    in_header = True
    current_motif = []

    for line in lines:
        if in_header:
            if line.startswith("MOTIF"):
                in_header = False
                current_motif = [line]
            else:
                header_lines.append(line)
        else:
            # in motif region
            if line.startswith("MOTIF"):
                # Start a new motif, push previous if any
                if current_motif:
                    motifs.append(current_motif)
                current_motif = [line]
            else:
                current_motif.append(line)

    # Append last motif if exists
    if current_motif:
        motifs.append(current_motif)

    if not motifs:
        raise ValueError("No motifs found (no lines starting with 'MOTIF').")

    return header_lines, motifs


def sanitize_motif_name(name: str) -> str:
    """
    Turn motif name into a safe filename fragment.
    """
    # Strip leading/trailing whitespace and replace spaces with underscores
    name = name.strip().replace(" ", "_")
    # Remove anything that is obviously bad for filenames
    name = re.sub(r"[^A-Za-z0-9._+-]", "_", name)
    if not name:
        name = "motif"
    return name


def motif_name_from_block(motif_lines: List[str]) -> str:
    """
    Extract motif name from the first 'MOTIF' line of the motif block.
    MEME v4 'MOTIF' line is typically:
       MOTIF <name> [optional stuff...]
    """
    first_line = motif_lines[0].strip()
    # Remove leading "MOTIF"
    if not first_line.startswith("MOTIF"):
        raise ValueError("Motif block does not start with 'MOTIF': " + first_line)
    parts = first_line.split()
    if len(parts) < 2:
        # No explicit name? Fallback to generic
        return "motif"
    return parts[1]


def make_single_meme_text(
    header_lines: List[str], motif_lines: List[str]
) -> str:
    """
    Construct the contents of a single-motif MEME file as text.
    Assumes header_lines already contain 'MEME version 4', ALPHABET, background, etc.
    """
    # Ensure exactly one blank line between header and motif (optional aesthetic)
    content = "".join(header_lines)
    if not content.endswith("\n\n"):
        if not content.endswith("\n"):
            content += "\n"
        content += "\n"
    content += "".join(motif_lines)
    return content


def write_motifs_to_tar(
    header_lines: List[str],
    motifs: List[List[str]],
    output_tar_gz: str,
    base_prefix: str = "",
) -> None:
    """
    Write each motif as its own MEME file into a tar.gz archive.
    Files are named: <base_prefix><sanitized_motif_name>.meme
    """
    # If base_prefix is given, make sure it ends with underscore or slash for clarity
    base_prefix = base_prefix or ""
    with tarfile.open(output_tar_gz, "w:gz") as tar:
        for idx, motif_lines in enumerate(motifs, start=1):
            motif_name = motif_name_from_block(motif_lines)
            safe_name = sanitize_motif_name(motif_name)
            filename = f"{base_prefix}{safe_name}.meme"

            meme_text = make_single_meme_text(header_lines, motif_lines)
            data = meme_text.encode("utf-8")

            info = tarfile.TarInfo(name=filename)
            info.size = len(data)
            tar.addfile(info, BytesIO(data))


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Split a multi-motif MEME v4 file into separate per-motif MEME files "
            "and package them into a .tar.gz archive."
        )
    )
    parser.add_argument(
        "input_meme",
        help="Path to input MEME file containing many motifs",
    )
    parser.add_argument(
        "output_tar_gz",
        help="Path to output tar.gz file that will contain separate MEME files",
    )
    parser.add_argument(
        "--prefix",
        default="",
        help="Optional filename prefix for each per-motif MEME file inside the tar (default: none)",
    )
    args = parser.parse_args()

    lines = read_file(args.input_meme)
    header_lines, motifs = split_header_and_motifs(lines)

    print(f"Found {len(motifs)} motifs.")
    write_motifs_to_tar(header_lines, motifs, args.output_tar_gz, base_prefix=args.prefix)
    print(f"Written all motifs into tar.gz: {args.output_tar_gz}")


if __name__ == "__main__":
    main()

