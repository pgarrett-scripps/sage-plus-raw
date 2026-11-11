//! Sequential spectrum reading from one owned file handle.
//!
//! [`RawFileReader::open_path`] parses the file's structure and then drops the
//! handle it read it through, so a caller that wants spectra opens the same
//! path a second time. [`SpectrumReader`] keeps that first handle instead and
//! reads every scan through it, which saves the second open and lets a
//! forward seek inside the read buffer keep the buffer rather than discard
//! it.
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
use std::io::{self, BufReader, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};

use crate::error::{Error, Result};
use crate::spectrum::{try_extract_spectrum_into, SpectrumHeader, SpectrumRecord};
use crate::RawFileReader;

/// Read buffer held by a [`SpectrumReader`]. Large enough that a run of
/// scans is usually read from one buffer fill, since consecutive scan
/// packets are consecutive in the file.
const BUFFER_BYTES: usize = 1 << 20;

/// A buffered reader whose absolute seeks are expressed as relative ones, so
/// a seek to a position already inside the buffer keeps the buffered bytes.
///
/// [`BufReader`]'s own [`Seek`] implementation discards the buffer on every
/// seek, including a seek forward by a few bytes, which is what per-scan
/// absolute seeking does to a sequential read. [`BufReader::seek_relative`]
/// keeps the buffer when the target is inside it, so this tracks the logical
/// position and converts.
struct ScanSource<R> {
    inner: BufReader<R>,
    pos: u64,
}

impl<R: Read + Seek> ScanSource<R> {
    fn new(inner: R) -> Self {
        Self {
            inner: BufReader::with_capacity(BUFFER_BYTES, inner),
            pos: 0,
        }
    }
}

impl<R: Read> Read for ScanSource<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let read = self.inner.read(buf)?;
        self.pos += read as u64;
        Ok(read)
    }
}

impl<R: Read + Seek> Seek for ScanSource<R> {
    fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
        let target = match pos {
            SeekFrom::Start(offset) => offset,
            SeekFrom::Current(delta) => self.pos.checked_add_signed(delta).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidInput, "seek before start of file")
            })?,
            SeekFrom::End(delta) => {
                // Only an end-relative seek has to ask the file how long it
                // is, which discards the buffer; `BinaryReader` caches the
                // length, so this runs once per reader.
                let end = self.inner.seek(SeekFrom::End(0))?;
                self.pos = end;
                end.checked_add_signed(delta).ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidInput, "seek before start of file")
                })?
            }
        };
        match i64::try_from(i128::from(target) - i128::from(self.pos)) {
            Ok(delta) => self.inner.seek_relative(delta)?,
            // Further than an i64 can express: an absolute seek, which the
            // buffer could not have covered anyway.
            Err(_) => {
                self.inner.seek(SeekFrom::Start(target))?;
            }
        }
        self.pos = target;
        Ok(target)
    }

    fn stream_position(&mut self) -> io::Result<u64> {
        Ok(self.pos)
    }
}

/// A `.raw` file opened once, read sequentially through the handle it was
/// opened with.
///
/// The decoded values are exactly those of [`crate::try_extract_spectrum`]
/// and [`crate::try_iter_spectra`] on the same file; this only changes how
/// many handles and buffer fills getting them costs. Scans may be read in
/// any order, but a file read front to back is what the buffer suits.
pub struct SpectrumReader {
    raw: RawFileReader,
    source: ScanSource<File>,
    path: PathBuf,
}

impl SpectrumReader {
    /// Open `path`, parse its structure, and keep the handle for reading
    /// scans.
    pub fn open_path(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let result = (|| {
            let mut source = ScanSource::new(File::open(path)?);
            let raw = RawFileReader::open(&mut source)?;
            Ok(Self {
                raw,
                source,
                path: path.to_path_buf(),
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
        try_extract_spectrum_into(
            &self.raw,
            &mut self.source,
            idx,
            include_profile,
            mz,
            intensity,
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    /// A `Cursor` that counts the seeks reaching it, to tell a buffered
    /// forward seek from one that discards the buffer.
    struct CountingCursor {
        inner: Cursor<Vec<u8>>,
        seeks: u32,
    }

    impl Read for CountingCursor {
        fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
            self.inner.read(buf)
        }
    }

    impl Seek for CountingCursor {
        fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
            self.seeks += 1;
            self.inner.seek(pos)
        }
    }

    fn source(len: usize) -> ScanSource<CountingCursor> {
        let bytes: Vec<u8> = (0..len).map(|i| i as u8).collect();
        ScanSource::new(CountingCursor {
            inner: Cursor::new(bytes),
            seeks: 0,
        })
    }

    fn byte_at(s: &mut ScanSource<CountingCursor>, offset: u64) -> u8 {
        s.seek(SeekFrom::Start(offset)).unwrap();
        let mut buf = [0u8; 1];
        s.read_exact(&mut buf).unwrap();
        buf[0]
    }

    #[test]
    fn every_seek_form_lands_where_the_bytes_say() {
        let mut s = source(4096);
        assert_eq!(byte_at(&mut s, 0), 0);
        assert_eq!(byte_at(&mut s, 300), 300u32 as u8);
        assert_eq!(s.stream_position().unwrap(), 301);
        s.seek(SeekFrom::Current(-101)).unwrap();
        assert_eq!(s.stream_position().unwrap(), 200);
        let mut buf = [0u8; 2];
        s.read_exact(&mut buf).unwrap();
        assert_eq!(buf, [200u32 as u8, 201u32 as u8]);
        assert_eq!(s.seek(SeekFrom::End(-1)).unwrap(), 4095);
        s.read_exact(&mut buf[..1]).unwrap();
        assert_eq!(buf[0], 4095u32 as u8);
        assert_eq!(s.stream_position().unwrap(), 4096);
    }

    #[test]
    fn a_forward_seek_inside_the_buffer_does_not_reach_the_file() {
        let mut s = source(4096);
        // One fill brings the whole file into the buffer.
        assert_eq!(byte_at(&mut s, 0), 0);
        let after_fill = s.inner.get_ref().seeks;
        // What a sequential scan loop does: absolute seek forward a little,
        // then read. Both stay inside the buffered bytes.
        for offset in [16u64, 64, 1000, 4000] {
            assert_eq!(byte_at(&mut s, offset), offset as u8);
        }
        assert_eq!(
            s.inner.get_ref().seeks,
            after_fill,
            "an in-buffer seek should not seek the file"
        );
        // Backwards past the buffer start does have to seek.
        s.seek(SeekFrom::End(0)).unwrap();
        assert!(s.inner.get_ref().seeks > after_fill);
    }

    #[test]
    fn a_seek_before_the_start_of_the_file_is_an_error() {
        let mut s = source(64);
        assert!(s.seek(SeekFrom::Current(-1)).is_err());
        assert!(s.seek(SeekFrom::End(-65)).is_err());
    }
}
