#!/usr/bin/env python
"""
set_singularity_bind_paths.py <config_filename> <samples_filename>
"""
import os
import pathlib
import sys
import yaml


def print_bind_paths(config_filename, samples_filename):
    bind_paths = resolve_additional_bind_paths(
        get_paths(config_filename, samples_filename)
    )
    print(",".join(bind_paths))


def get_paths(config_filename, samples_filename):
    paths = { config_filename, samples_filename }
    with open(config_filename, "r") as config_file:
        config = yaml.safe_load(config_file)
    conf_keys = [
        "workdir",
        "scriptsdir",
        "samplemanifest",
        "contrasts",
        "preparsedDir",
        "adapters",
        "ccbr_tools_path"
    ]
    paths.update({ config[key] for key in conf_keys if key in config })
    paths.update(config['reference'][config['genome']].values())
    paths.update(config['spikein_reference'][config['spikein_genome']].values())

    with open(samples_filename, "r") as samples_file:
        next(samples_file)  # skip header
        for line in samples_file:
            line_spl = line.split("\t")
            if len(line_spl) > 1:
                paths.add(line_spl[1])
            if len(line_spl) > 2:
                paths.add(line_spl[2])

    return paths


def resolve_additional_bind_paths(search_paths):
    """Adapted from RENEE

    Finds additional singularity bind paths from a list of random paths. Paths are
    indexed with a compostite key containing the first two directories of an absolute
    file path to avoid issues related to shared names across the /gpfs shared network
    filesystem. For each indexed list of file paths, a common path is found. Assumes
    that the paths provided are absolute paths, the renee build sub command creates
    resource file index with absolute filenames.
    @param search_paths list[<str>]:
        List of absolute file paths to find common bind paths from
    @return common_paths list[<str>]:
        Returns a list of common shared file paths to create additional singularity bind paths
    """
    common_paths = []
    indexed_paths = {}

    for ref in search_paths:
        # Skip over resources with remote URI and
        # skip over strings that are not file PATHS as
        # RENEE build creates absolute resource PATHS
        if (
            ref.lower().startswith("sftp://")
            or ref.lower().startswith("s3://")
            or ref.lower().startswith("gs://")
            or not ref.lower().startswith(os.sep)
        ):
            continue

        # Break up path into directory tokens
        for r in [
            ref,
            str(pathlib.Path(ref).resolve()),
        ]:  # taking care of paths which are symlinks!
            path_list = os.path.abspath(r).split(os.sep)

            try:  # Create composite index from first two directories
                # Avoids issues created by shared /gpfs/ PATHS
                index = path_list[1:3]
                index = tuple(index)
            except IndexError:
                index = path_list[1]  # ref startswith /
            if index not in indexed_paths:
                indexed_paths[index] = []
            # Create an INDEX to find common PATHS for each root child directory
            # like /scratch or /data. This prevents issues when trying to find the
            # common path between these two different directories (resolves to /)
            indexed_paths[index].append(str(os.sep).join(path_list))

    for index, paths in indexed_paths.items():
        # Find common paths for each path index
        common_paths.append(os.path.dirname(os.path.commonprefix(paths)))

    return sorted(set(common_paths))


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__)
        raise Exception("Please provide two arguments")
    config_filename = sys.argv[1]
    samples_filename = sys.argv[2]
    print_bind_paths(config_filename, samples_filename)
