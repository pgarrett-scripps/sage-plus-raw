# Credits

## Upstream OpenTFRaw

This independently maintained fork derives from [OpenTFRaw](https://github.com/Sigilweaver/OpenTFRaw), originally authored by Nathan Riley and contributors. Original copyright 2026 Sigilweaver Holdings LLC and Apache-2.0 licensing are retained. See NOTICE, LICENSE and CONTRIBUTORS.md. The fork does not own the upstream opentfraw registry identity or website.

Qualification also uses publicly distributed ProteoWizard RAW and reference mzML fixtures at the revision recorded in scripts/public-fixtures.json. Those files are downloaded and hash-checked during tests, not redistributed in the Rust package.

## Prior art

### unfinnigan

Gene Selkov, 2010-2012. Perl and Python reverse-engineering of the Thermo RAW binary format.
The most thorough prior independent analysis of the format, covering versions 57, 62, 63, 64,
and 66. Field names and layout notes from unfinnigan were cross-referenced when validating
field offsets.

Source: https://github.com/prvst/unfinnigan

## Format specifications

The embedded method file is read as a Compound File Binary container, following Microsoft's
open specification [MS-CFB](https://learn.microsoft.com/openspecs/windows_protocols/ms-cfb/)
("Compound File Binary File Format"). Stream reads were checked during development against
[olefile](https://github.com/decalage2/olefile) 0.47 (Philippe Lagadec, BSD-2-Clause), an
independent open-source reader of that format. It is not a dependency.

## Rust dependencies

- [roxmltree](https://github.com/RazrFalcon/roxmltree) -- read-only parsing of
  embedded method XML (Yevhenii Reizner, MIT/Apache-2.0).

- [thiserror](https://github.com/dtolnay/thiserror) -- derive macro for Error impls (David Tolnay, MIT/Apache-2.0)
