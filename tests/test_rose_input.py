import contextlib
import importlib.util
import io
import tempfile
import unittest
from pathlib import Path
from unittest import mock

SCRIPT = Path(__file__).resolve().parents[1] / "workflow/scripts/_prep_rose_input.py"
SPEC = importlib.util.spec_from_file_location("prep_rose_input", SCRIPT)
rose_input = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(rose_input)


class TestRoseInput(unittest.TestCase):
    def test_numeric_scores_do_not_require_math_isfinite(self):
        with mock.patch.object(rose_input.math, "isfinite", side_effect=AttributeError):
            self.assertEqual(rose_input.parse_number("1.25"), 1.25)
            self.assertEqual(rose_input.parse_number("-2"), -2.0)

    def test_invalid_and_nonfinite_scores_warn_and_use_default(self):
        for value in ("not-a-number", None, "NaN", "inf", "-inf"):
            with (
                self.subTest(value=value),
                contextlib.redirect_stderr(io.StringIO()) as err,
            ):
                self.assertEqual(rose_input.parse_number(value, default=7.0), 7.0)
                self.assertIn("peak score", err.getvalue())

    def test_peak_bed_and_gff_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)
            peaks = path / "peaks.bed"
            peaks.write_text("chr1\t5\t10\tpeak\t2.5\t+\n")
            rows = rose_input.load_peaks(str(peaks), "bed", "sample")
            bed = path / "output.bed"
            gff = path / "output.gff"
            rose_input.write_bed6(rows, str(bed))
            rose_input.convert_stitched_bed_to_rose_gff(str(bed), str(gff))
            self.assertEqual(
                bed.read_text(), "chr1\t5\t10\tsample_000001\t2.500000\t+\n"
            )
            self.assertEqual(
                gff.read_text(),
                "chr1\tsample_000001\t.\t6\t10\t.\t.\t.\t.\t.\t.\tsample_000001\n",
            )


if __name__ == "__main__":
    unittest.main()
