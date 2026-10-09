//! Native Thermo RAW reading with strict spectrum extraction.
//!
//! Cargo package: `opentfraw-reader`. Rust library name: `opentfraw`.
//! This is an independent OpenTFRaw fork. The launch candidate is not yet
//! published to crates.io. See the packaged README and REVIEW for support scope.
//!
//! ```no_run
//! use opentfraw::{try_iter_spectra, RawFileReader};
//! use std::{fs::File, io::BufReader};
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let raw = RawFileReader::open_path("acquisition.raw")?;
//! let mut source = BufReader::new(File::open("acquisition.raw")?);
//! for spectrum in try_iter_spectra(&raw, &mut source, false) {
//!     let spectrum = spectrum?;
//!     println!("scan {}: {} peaks", spectrum.scan_number, spectrum.mz.len());
//! }
//! # Ok(())
//! # }
//! ```
//!
//! Strict extraction propagates unreadable arrays. Unknown scientific metadata
//! remains unknown. Convenience mzML export rejects unsupported precursor modes
//! and unqualified activation. Discard partial output if a writer returns an error.
//! The legacy extraction helpers discard errors and exist only for compatibility.

#![cfg_attr(not(test), warn(clippy::unwrap_used, clippy::expect_used))]
mod binary_reader;
mod error;
mod event_stream;
mod metadata;
mod reader;
mod spectrum;
mod types;

pub(crate) mod bytes;
#[cfg(test)]
mod test_util;

pub mod audit_tag;
pub mod device;
pub mod error_log;
pub mod extra;
pub mod generic_data;
pub mod header;
pub mod mzml;
pub mod raw_file_info;
pub mod run_header;
pub mod sample_info;
pub mod scan_data;
pub mod scan_event;
pub mod scan_filter;
pub mod scan_format;
pub mod scan_index;
pub mod seq_row;

pub use device::{DetectedInstrument, DeviceFamily};
pub use error::{Error, ErrorContext, ErrorKind, Operation, Result};
pub use extra::ExtraFields;
pub use metadata::{ScanParams, StatusLogEntry};
#[allow(deprecated)]
pub use mzml::{
    extract_spectrum, iter_spectra, scan_metadata, try_extract_spectrum,
    try_extract_spectrum_arrays, try_iter_spectra, write_indexed_mzml, write_mzml, PrecursorInfo,
    ScanMetadata, SpectrumArrays, SpectrumRecord,
};
pub use reader::{ControllerInfo, ControllerType, RawFileReader};
pub use scan_format::ScanDataFormat;
pub use types::*;
