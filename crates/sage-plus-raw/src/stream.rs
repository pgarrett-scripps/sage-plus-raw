//! Sequential spectrum reading from one owned file handle.
//!
//! [`RawFileReader::open_path`] parses the file's structure and then drops
//! the handle it read it through, so a caller that wants spectra opens the
//! same path a second time and carries both. [`SpectrumReader`] keeps that
//! first handle and reads every scan through it.
//!
//! This is a convenience, not a faster decoder. Measured over 67,599 scans,
//! the second open is one syscall, and a scan loop spends its time on the
//! roughly two reads per peak that decoding a packet does, around 2,400 per
//! scan. Per-scan seeking was cut where cutting it pays, in
//! [`RawFileReader::file_len`], which this and the existing path both read
//! through; interposing anything on the reads themselves cost more than it
//! saved, so this holds a plain [`BufReader`].
//!
//! ```no_run
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! for spectrum in sage_plus_raw::iter_spectra_from_path("acquisition.raw", false)? {
//!     let spectrum = spectrum?;
//!     println!("scan {}: {} peaks", spectrum.scan_number, spectrum.mz.len());
//! }
//! # Ok(())
//! # }
//! ```

use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};

use crate::error::{Error, Result};
use crate::spectrum::{SpectrumHeader, SpectrumRecord};
use crate::RawFileReader;

/// A `.raw` file opened once, read through the handle it was opened with.
///
/// The decoded values are exactly those of [`crate::try_extract_spectrum`]
/// and [`crate::try_iter_spectra`] on the same file.
pub struct SpectrumReader {
    raw: RawFileReader,
    source: BufReader<File>,
    path: PathBuf,
    /// Reused across scans: each scan's peak block is read into it in one
    /// go and parsed from the bytes.
    scratch: Vec<u8>,
}

impl SpectrumReader {
    /// Open `path`, parse its structure, and keep the handle for reading
    /// scans.
    pub fn open_path(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let result = (|| {
            let mut source = BufReader::new(File::open(path)?);
            let raw = RawFileReader::open(&mut source)?;
            Ok(Self {
                raw,
                source,
                path: path.to_path_buf(),
                scratch: Vec::new(),
            })
        })();
        result.map_err(|error: Error| {
            let mut context = crate::ErrorContext::new(crate::Operation::Open);
            context.path = Some(path.to_path_buf());
            error.with_context(context)
        })
    }

    /// The parsed file structure: scan index, scan events, run header and
    /// instrument metadata.
    pub fn raw(&self) -> &RawFileReader {
        &self.raw
    }

    /// The path this reader was opened from.
    pub fn path(&self) -> &Path {
        &self.path
    }

    /// Number of scans in the file.
    pub fn num_scans(&self) -> u32 {
        self.raw.num_scans
    }

    /// Decode the scan at `idx` into a new record, as
    /// [`crate::try_extract_spectrum`] does.
    pub fn read_spectrum(&mut self, idx: u32, include_profile: bool) -> Result<SpectrumRecord> {
        let mut mz = Vec::new();
        let mut intensity = Vec::new();
        let header = self.read_spectrum_into(idx, include_profile, &mut mz, &mut intensity)?;
        Ok(header.into_record(mz, intensity))
    }

    /// Decode the scan at `idx` into caller-owned buffers, as
    /// [`crate::try_extract_spectrum_into`] does: `mz` and `intensity` are
    /// cleared and refilled, keeping their capacity, and are left empty on an
    /// error.
    pub fn read_spectrum_into(
        &mut self,
        idx: u32,
        include_profile: bool,
        mz: &mut Vec<f64>,
        intensity: &mut Vec<f32>,
    ) -> Result<SpectrumHeader> {
        crate::spectrum::try_extract_spectrum_into_with(
            &self.raw,
            &mut self.source,
            idx,
            include_profile,
            mz,
            intensity,
            &mut self.scratch,
        )
    }

    /// Iterate every scan without dropping failures, as
    /// [`crate::try_iter_spectra`] does. Each item corresponds to exactly one
    /// scan index.
    pub fn iter_spectra(
        &mut self,
        include_profile: bool,
    ) -> impl Iterator<Item = Result<SpectrumRecord>> + '_ {
        (0..self.raw.num_scans).map(move |idx| self.read_spectrum(idx, include_profile))
    }

    /// Turn this reader into an owning iterator over every scan, for a
    /// caller that wants the iterator to outlive the local binding.
    pub fn into_spectra(self, include_profile: bool) -> SpectrumIter {
        SpectrumIter {
            reader: self,
            include_profile,
            next_idx: 0,
        }
    }
}

/// Owning iterator over every scan of one file, from
/// [`SpectrumReader::into_spectra`] or [`iter_spectra_from_path`].
///
/// Like [`crate::try_iter_spectra`], it yields one item per scan index and
/// does not stop at a failed scan; collect into `Result<Vec<_>>` to stop on
/// the first error.
pub struct SpectrumIter {
    reader: SpectrumReader,
    include_profile: bool,
    next_idx: u32,
}

impl SpectrumIter {
    /// The reader underneath, for the file's metadata.
    pub fn reader(&self) -> &SpectrumReader {
        &self.reader
    }

    /// Give the reader back, to read scans by index or into buffers.
    pub fn into_reader(self) -> SpectrumReader {
        self.reader
    }
}

impl Iterator for SpectrumIter {
    type Item = Result<SpectrumRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.next_idx >= self.reader.raw.num_scans {
            return None;
        }
        let idx = self.next_idx;
        self.next_idx += 1;
        Some(self.reader.read_spectrum(idx, self.include_profile))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let left = self.reader.raw.num_scans.saturating_sub(self.next_idx) as usize;
        (left, Some(left))
    }
}

impl ExactSizeIterator for SpectrumIter {}

/// Open `path` and iterate every scan from that one handle.
///
/// Equivalent to [`RawFileReader::open_path`] plus
/// [`crate::try_iter_spectra`] over a second handle on the same path, without
/// the second open.
pub fn iter_spectra_from_path(
    path: impl AsRef<Path>,
    include_profile: bool,
) -> Result<SpectrumIter> {
    Ok(SpectrumReader::open_path(path)?.into_spectra(include_profile))
}
