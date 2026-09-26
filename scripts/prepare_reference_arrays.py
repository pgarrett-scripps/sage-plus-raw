"""Normalize verified independent mzML arrays for dependency-free Rust tests.

Run fetch_public_fixtures.py first. Only SHA-256-pinned references are accepted.
The binary test format is scan count, then index, MS level, RT in minutes, count,
and paired f64 m/z and intensity values. Integers and floats are little endian.
"""
import base64
import hashlib
import json
from pathlib import Path
import struct
import xml.etree.ElementTree as ET
import zlib

root = Path(__file__).resolve().parent.parent
manifest = json.loads((root / "scripts/public-fixtures.json").read_text())
ns = {"m": "http://psi.hupo.org/ms/mzml"}
for name, digest in manifest["files"].items():
    if not name.endswith(".mzML"):
        continue
    path = root / "corpus/public" / name
    content = path.read_bytes()
    if hashlib.sha256(content).hexdigest() != digest:
        raise ValueError(f"Reference digest mismatch: {name}")
    spectra = ET.fromstring(content).findall(".//m:spectrum", ns)
    output = bytearray(struct.pack("<I", len(spectra)))
    for spectrum in spectra:
        cv = {x.attrib["accession"]: x.attrib for x in spectrum.findall(".//m:cvParam", ns)}
        level = int(cv["MS:1000511"]["value"])
        time = cv["MS:1000016"]
        rt = float(time["value"])
        if time["unitAccession"] == "UO:0000010":
            rt /= 60
        elif time["unitAccession"] != "UO:0000031":
            raise ValueError("Unknown retention time unit")
        arrays = {}
        for array in spectrum.findall("m:binaryDataArrayList/m:binaryDataArray", ns):
            terms = {x.attrib["accession"] for x in array.findall("m:cvParam", ns)}
            key = "mz" if "MS:1000514" in terms else "intensity" if "MS:1000515" in terms else None
            if key is None:
                continue
            data = base64.b64decode(array.findtext("m:binary", namespaces=ns))
            if "MS:1000574" in terms:
                data = zlib.decompress(data)
            elif "MS:1000576" not in terms:
                raise ValueError("Unknown compression")
            dtype = "d" if "MS:1000523" in terms else "f" if "MS:1000521" in terms else None
            if dtype is None:
                raise ValueError("Unknown numeric representation")
            if dtype != "f":
                raise ValueError("Qualification tolerances require 32-bit reference arrays")
            arrays[key] = [x[0] for x in struct.iter_unpack("<" + dtype, data)]
        count = int(spectrum.attrib["defaultArrayLength"])
        if len(arrays["mz"]) != count or len(arrays["intensity"]) != count:
            raise ValueError("Reference array length mismatch")
        output += struct.pack("<IIdI", int(spectrum.attrib["index"]), level, rt, count)
        for mz, intensity in zip(arrays["mz"], arrays["intensity"]):
            output += struct.pack("<dd", mz, intensity)
    path.with_suffix(".arrays").write_bytes(output)
