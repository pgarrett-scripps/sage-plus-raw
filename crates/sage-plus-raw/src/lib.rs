//! Native Thermo RAW reading with strict spectrum extraction, for Sage Plus.
//!
//! Fork of Sigilweaver OpenTFRaw by Nathan Riley and contributors; not
//! affiliated with Sigilweaver or Thermo Fisher Scientific.
//!
//! ```no_run
//! use sage_plus_raw::{try_iter_spectra, RawFileReader};
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
//! Strict extraction propagates unreadable scans. Invalid single peaks are
//! dropped and counted in `dropped_peaks`. Unknown scientific metadata
//! remains unknown.

#![cfg_attr(not(test), warn(clippy::unwrap_used, clippy::expect_used))]
mod binary_reader;
mod error;
mod event_stream;
mod instrument_method;
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
pub mod generic_data;
pub mod header;
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
pub use metadata::{ScanParams, StatusLogEntry};
pub use reader::{ControllerInfo, ControllerType, RawFileReader};
pub use scan_format::ScanDataFormat;
pub use spectrum::{
    scan_metadata, try_extract_spectrum, try_extract_spectrum_arrays, try_iter_spectra,
    PrecursorInfo, ScanMetadata, SpectrumArrays, SpectrumRecord,
};
pub use types::*;
