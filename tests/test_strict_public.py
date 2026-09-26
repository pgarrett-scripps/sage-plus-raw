"""Runtime checks for the strict Python boundary using pinned public fixtures."""
from pathlib import Path

import numpy as np
import pytest
import opentfraw

ROOT = Path(__file__).resolve().parent.parent / "corpus/public"


def test_centroid_scan_mode_describes_returned_arrays():
    raw = opentfraw.RawFile(str(ROOT / "BSA-FT-HCD.raw"))
    scan = raw.scan(raw.first_scan)
    assert scan["scan_mode"] == "centroid"
    assert scan["intensity"].dtype == np.float32
    mz, intensity = raw.profile(raw.first_scan)
    assert len(mz) == len(intensity) == 509
    assert intensity.dtype == np.float64
    assert np.isfinite(mz).all() and np.isfinite(intensity).all()


def test_errors_have_distinct_python_categories(tmp_path):
    with pytest.raises(OSError):
        opentfraw.RawFile(str(tmp_path / "missing.raw"))
    truncated = tmp_path / "truncated.raw"
    truncated.write_bytes(b"")
    with pytest.raises(ValueError, match="Open"):
        opentfraw.RawFile(str(truncated))
    raw = opentfraw.RawFile(str(ROOT / "BSA-FT-HCD.raw"))
    with pytest.raises(IndexError):
        raw.peaks(raw.last_scan + 1)
    multiplexed = opentfraw.RawFile(str(ROOT / "FT-HCD-MSX.raw"))
    with pytest.raises(NotImplementedError, match="multiple precursor"):
        multiplexed.to_mzml(str(tmp_path / "unsupported.mzML"))


def test_export_is_complete_on_return(tmp_path):
    raw = opentfraw.RawFile(str(ROOT / "BSA-FT-HCD.raw"))
    output = tmp_path / "qualified.mzML"
    raw.to_mzml(str(output))
    assert output.read_text().rstrip().endswith("</mzML>")
