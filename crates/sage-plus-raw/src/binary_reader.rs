//! Bounded binary I/O used by RAW structure decoders.
use std::io::{Read, Seek, SeekFrom};

use crate::error::{Error, Result};

/// Low-level binary reading helpers.
pub(crate) struct BinaryReader<R> {
    inner: R,
    pos: u64,
    /// True offset in the underlying stream that `pos == 0` means, when the
    /// caller knew it. `None` means it has to be asked for. Kept separate
    /// from `pos` so reported error offsets stay relative to where this
    /// reader started, as they have always been.
    base: Option<u64>,
    /// Cached total stream length, populated on first call to [`Self::length`]
    /// or [`Self::remaining`]. The underlying stream is never written to
    /// while parsing, so this is safe to cache for the reader's lifetime and
    /// avoids repeated end-of-stream seeks in allocation-size checks that can
    /// run once per record in a large loop.
    total_len: Option<u64>,
}

impl<R: Read + Seek> BinaryReader<R> {
    pub fn new(inner: R) -> Self {
        Self {
            inner,
            pos: 0,
            base: None,
            total_len: None,
        }
    }

    /// Like [`Self::new`], but told where in the stream it starts and how
    /// long the stream is, so bounded reads can check what is left without
    /// asking the stream anything.
    ///
    /// [`Self::new`] has to ask: it seeks to the end and back for the
    /// length, then asks the position again for every bounded read. A reader
    /// built per scan pays all of that per scan.
    pub fn at(inner: R, base: u64, total_len: u64) -> Self {
        Self {
            inner,
            pos: 0,
            base: Some(base),
            total_len: Some(total_len),
        }
    }

    pub fn into_inner(self) -> R {
        self.inner
    }

    #[allow(dead_code)]
    pub(crate) fn position(&self) -> u64 {
        self.pos
    }

    pub fn seek_to(&mut self, offset: u64) -> Result<()> {
        self.inner.seek(SeekFrom::Start(offset))?;
        self.pos = offset;
        // `pos` now counts from the start of the stream, not from where this
        // reader began, so a recorded base no longer describes it.
        self.base = None;
        Ok(())
    }

    pub fn read_bytes(&mut self, n: usize) -> Result<Vec<u8>> {
        // Guard against a crafted/corrupt file declaring an implausible
        // length-prefix: verify the bytes could actually exist in the
        // remaining input before allocating a buffer for them.
        self.check_count(n as u64, 1)?;
        let mut buf = vec![0u8; n];
        self.inner.read_exact(&mut buf).map_err(|e| {
            if e.kind() == std::io::ErrorKind::UnexpectedEof {
                Error::UnexpectedEof {
                    offset: self.pos,
                    needed: n,
                }
            } else {
                Error::Io(e)
            }
        })?;
        self.pos += n as u64;
        Ok(buf)
    }

    pub fn read_bytes_into(&mut self, buf: &mut [u8]) -> Result<()> {
        let n = buf.len();
        self.inner.read_exact(buf).map_err(|e| {
            if e.kind() == std::io::ErrorKind::UnexpectedEof {
                Error::UnexpectedEof {
                    offset: self.pos,
                    needed: n,
                }
            } else {
                Error::Io(e)
            }
        })?;
        self.pos += n as u64;
        Ok(())
    }

    pub fn skip(&mut self, n: usize) -> Result<()> {
        self.inner.seek(SeekFrom::Current(n as i64))?;
        self.pos += n as u64;
        Ok(())
    }

    pub fn length(&mut self) -> Result<u64> {
        if let Some(len) = self.total_len {
            return Ok(len);
        }
        // Query and restore the *true* underlying position (not `self.pos`):
        // some callers construct a fresh `BinaryReader` over a source that
        // was already seeked to a nonzero offset, so `self.pos` (which starts
        // at 0 for every new reader) would not agree with where the
        // underlying stream actually is.
        let cur = self.inner.stream_position()?;
        let end = self.inner.seek(SeekFrom::End(0))?;
        self.inner.seek(SeekFrom::Start(cur))?;
        self.total_len = Some(end);
        Ok(end)
    }

    /// Bytes remaining between the stream's true current position and its end.
    pub(crate) fn remaining(&mut self) -> Result<u64> {
        if let (Some(base), Some(end)) = (self.base, self.total_len) {
            return Ok(end.saturating_sub(base.saturating_add(self.pos)));
        }
        let cur = self.inner.stream_position()?;
        let end = self.length()?;
        Ok(end.saturating_sub(cur))
    }

    /// Verify that `count` items of at least `min_item_bytes` bytes each
    /// could plausibly still be read from the stream, before any allocation
    /// proportional to `count` is attempted.
    ///
    /// Crafted or corrupt RAW files can declare implausible counts (a scan
    /// count derived from a bogus header field, a peak count with no
    /// relation to the file's actual size, etc). Without this check, code
    /// like `Vec::with_capacity(count)` would allocate memory proportional
    /// to an attacker-controlled value rather than to what the input could
    /// actually contain - a memory-exhaustion vector reachable before a
    /// single byte of the declared payload is read.
    pub(crate) fn check_count(&mut self, count: u64, min_item_bytes: u64) -> Result<()> {
        let remaining = self.remaining()?;
        let needed = count.saturating_mul(min_item_bytes.max(1));
        if needed > remaining {
            return Err(Error::AllocationTooLarge {
                offset: self.pos,
                requested: needed,
                available: remaining,
            });
        }
        Ok(())
    }

    /// Read `count` fixed-size items into `buf` in one go, for a stream of
    /// identical records such as a scan's peak list or profile signal, and
    /// return the bytes read.
    ///
    /// `buf` is a caller-owned scratch buffer: it grows to fit the largest
    /// run it is asked for and is then written over, never cleared, so a
    /// caller reading scan after scan neither allocates nor zeroes per scan.
    /// Only `buf[..count * item_bytes]` is this run's data, which is what
    /// comes back. The bound is the same [`Self::check_count`] a
    /// field-at-a-time loop applies, and a stream that ends early reports
    /// the item that could not be completed, at the offset that loop would
    /// have reached, so the error a truncated file produces does not change.
    pub(crate) fn read_items_into<'a>(
        &mut self,
        buf: &'a mut Vec<u8>,
        count: usize,
        item_bytes: usize,
    ) -> Result<&'a [u8]> {
        self.check_count(count as u64, item_bytes as u64)?;
        let total = count.saturating_mul(item_bytes);
        if buf.len() < total {
            buf.resize(total, 0);
        }
        let mut filled = 0;
        while filled < total {
            match self.inner.read(&mut buf[filled..total]) {
                Ok(0) => {
                    // Stop where the loop would have: at the start of the
                    // first item it could not read whole.
                    let whole = filled - filled % item_bytes.max(1);
                    self.pos += whole as u64;
                    return Err(Error::UnexpectedEof {
                        offset: self.pos,
                        needed: item_bytes,
                    });
                }
                Ok(read) => filled += read,
                Err(error) if error.kind() == std::io::ErrorKind::Interrupted => continue,
                Err(error) => return Err(Error::Io(error)),
            }
        }
        self.pos += total as u64;
        Ok(&buf[..total])
    }

    pub fn read_u8(&mut self) -> Result<u8> {
        let mut buf = [0u8; 1];
        self.read_bytes_into(&mut buf)?;
        Ok(buf[0])
    }

    pub fn read_u16(&mut self) -> Result<u16> {
        let mut buf = [0u8; 2];
        self.read_bytes_into(&mut buf)?;
        Ok(u16::from_le_bytes(buf))
    }

    pub fn read_i16(&mut self) -> Result<i16> {
        let mut buf = [0u8; 2];
        self.read_bytes_into(&mut buf)?;
        Ok(i16::from_le_bytes(buf))
    }

    pub fn read_u32(&mut self) -> Result<u32> {
        let mut buf = [0u8; 4];
        self.read_bytes_into(&mut buf)?;
        Ok(u32::from_le_bytes(buf))
    }

    pub fn read_i32(&mut self) -> Result<i32> {
        let mut buf = [0u8; 4];
        self.read_bytes_into(&mut buf)?;
        Ok(i32::from_le_bytes(buf))
    }

    pub fn read_u64(&mut self) -> Result<u64> {
        let mut buf = [0u8; 8];
        self.read_bytes_into(&mut buf)?;
        Ok(u64::from_le_bytes(buf))
    }

    pub fn read_f32(&mut self) -> Result<f32> {
        let mut buf = [0u8; 4];
        self.read_bytes_into(&mut buf)?;
        Ok(f32::from_le_bytes(buf))
    }

    pub fn read_f64(&mut self) -> Result<f64> {
        let mut buf = [0u8; 8];
        self.read_bytes_into(&mut buf)?;
        Ok(f64::from_le_bytes(buf))
    }

    pub fn read_i8(&mut self) -> Result<i8> {
        let mut buf = [0u8; 1];
        self.read_bytes_into(&mut buf)?;
        Ok(buf[0] as i8)
    }

    /// Read a fixed-width UTF-16-LE string of `byte_len` bytes, stripping null padding.
    pub fn read_utf16_fixed(&mut self, byte_len: usize) -> Result<String> {
        let pos = self.pos;
        let raw = self.read_bytes(byte_len)?;
        if byte_len % 2 != 0 {
            return Err(Error::InvalidUtf16(pos));
        }
        let units: Vec<u16> = raw
            .chunks_exact(2)
            .map(|c| u16::from_le_bytes([c[0], c[1]]))
            .collect();
        // Find null terminator
        let end = units.iter().position(|&u| u == 0).unwrap_or(units.len());
        String::from_utf16(&units[..end]).map_err(|_| Error::InvalidUtf16(pos))
    }

    /// Like [`Self::read_utf16_fixed`], but replaces invalid UTF-16 (e.g.
    /// unpaired surrogates) with U+FFFD instead of failing. Used for
    /// free-text values such as trailer strings, where one bad character
    /// should not make the whole file unreadable.
    pub fn read_utf16_fixed_lossy(&mut self, byte_len: usize) -> Result<String> {
        let pos = self.pos;
        let raw = self.read_bytes(byte_len)?;
        if byte_len % 2 != 0 {
            return Err(Error::InvalidUtf16(pos));
        }
        let units: Vec<u16> = raw
            .chunks_exact(2)
            .map(|c| u16::from_le_bytes([c[0], c[1]]))
            .collect();
        let end = units.iter().position(|&u| u == 0).unwrap_or(units.len());
        Ok(String::from_utf16_lossy(&units[..end]))
    }

    /// Read a PascalStringWin32: UInt32 char count, then that many UTF-16-LE code units.
    pub fn read_pascal_string(&mut self) -> Result<String> {
        let pos = self.pos;
        let char_count = self.read_u32()? as usize;
        if char_count == 0 {
            return Ok(String::new());
        }
        let byte_len = char_count.checked_mul(2).ok_or(Error::InvalidUtf16(pos))?;
        let raw = self.read_bytes(byte_len)?;
        let units: Vec<u16> = raw
            .chunks_exact(2)
            .map(|c| u16::from_le_bytes([c[0], c[1]]))
            .collect();
        // Strip trailing nulls
        let end = units.iter().position(|&u| u == 0).unwrap_or(units.len());
        String::from_utf16(&units[..end]).map_err(|_| Error::InvalidUtf16(pos))
    }

    /// Read a Windows FILETIME and return Unix timestamp as f64 seconds.
    pub fn read_windows_filetime(&mut self) -> Result<f64> {
        let ft = self.read_u64()?;
        if ft == 0 {
            return Ok(0.0);
        }
        Ok((ft as f64 / 10_000_000.0) - 11_644_473_600.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn read_bytes_returns_requested_slice() {
        let mut r = BinaryReader::new(Cursor::new(vec![1, 2, 3, 4, 5]));
        assert_eq!(r.read_bytes(3).unwrap(), vec![1, 2, 3]);
        assert_eq!(r.read_bytes(2).unwrap(), vec![4, 5]);
    }

    #[test]
    fn read_bytes_past_end_is_eof_not_a_huge_allocation() {
        let mut r = BinaryReader::new(Cursor::new(vec![1u8, 2, 3]));
        let err = r.read_bytes(10).unwrap_err();
        assert!(matches!(err, Error::AllocationTooLarge { .. }));
    }

    #[test]
    fn read_bytes_rejects_implausible_declared_length() {
        // A corrupt/crafted file declaring an enormous length must be
        // rejected before any allocation is attempted, not after.
        let mut r = BinaryReader::new(Cursor::new(vec![0u8; 16]));
        let err = r.read_bytes(1_usize << 40).unwrap_err();
        match err {
            Error::AllocationTooLarge {
                requested,
                available,
                ..
            } => {
                assert_eq!(requested, 1u64 << 40);
                assert_eq!(available, 16);
            }
            other => panic!("expected AllocationTooLarge, got {other:?}"),
        }
    }

    #[test]
    fn check_count_allows_exact_fit() {
        let mut r = BinaryReader::new(Cursor::new(vec![0u8; 16]));
        assert!(r.check_count(2, 8).is_ok());
    }

    #[test]
    fn check_count_rejects_when_short_by_one_byte() {
        let mut r = BinaryReader::new(Cursor::new(vec![0u8; 15]));
        assert!(r.check_count(2, 8).is_err());
    }

    #[test]
    fn check_count_accounts_for_current_position() {
        let mut r = BinaryReader::new(Cursor::new(vec![0u8; 16]));
        r.skip(10).unwrap();
        // Only 6 bytes remain; asking for 1 item of 8 bytes must fail.
        assert!(r.check_count(1, 8).is_err());
        assert!(r.check_count(1, 6).is_ok());
    }

    #[test]
    fn check_count_saturates_instead_of_overflowing() {
        let mut r = BinaryReader::new(Cursor::new(vec![0u8; 4]));
        // count * min_item_bytes would overflow a u64 multiply; this must
        // not panic and must still be reported as too large.
        assert!(r.check_count(u64::MAX, u64::MAX).is_err());
    }

    #[test]
    fn length_is_cached_across_calls() {
        let mut r = BinaryReader::new(Cursor::new(vec![0u8; 100]));
        assert_eq!(r.length().unwrap(), 100);
        r.skip(50).unwrap();
        // Still reports the whole-stream length, not affected by position.
        assert_eq!(r.length().unwrap(), 100);
    }

    #[test]
    fn pascal_string_with_implausible_char_count_is_rejected_before_allocating() {
        let mut bytes = u32::MAX.to_le_bytes().to_vec(); // char_count
        bytes.extend_from_slice(&[0u8; 8]); // far short of the declared length
        let mut r = BinaryReader::new(Cursor::new(bytes));
        assert!(r.read_pascal_string().is_err());
    }

    #[test]
    fn utf16_fixed_strips_null_padding() {
        let mut bytes = Vec::new();
        for u in "hi".encode_utf16() {
            bytes.extend_from_slice(&u.to_le_bytes());
        }
        bytes.extend_from_slice(&[0, 0, 0, 0]); // null padding
        let mut r = BinaryReader::new(Cursor::new(bytes));
        assert_eq!(r.read_utf16_fixed(8).unwrap(), "hi");
    }

    #[test]
    fn utf16_fixed_stays_strict_about_unpaired_surrogates() {
        // Only schema-driven trailer values decode lossily
        // (`read_utf16_fixed_lossy`, used by `GenericRecord::read`). Fixed
        // structural strings keep rejecting invalid UTF-16, so a corrupt
        // structural field is still an error rather than a silent U+FFFD.
        let mut bytes = Vec::new();
        for u in [0x0041u16, 0xD800, 0x0042, 0x0000] {
            bytes.extend_from_slice(&u.to_le_bytes());
        }
        let mut r = BinaryReader::new(Cursor::new(bytes.clone()));
        assert!(matches!(r.read_utf16_fixed(8), Err(Error::InvalidUtf16(0))));
        let mut r = BinaryReader::new(Cursor::new(bytes));
        assert_eq!(r.read_utf16_fixed_lossy(8).unwrap(), "A\u{FFFD}B");
    }

    /// A `Cursor` that counts the seeks reaching it.
    struct CountingCursor {
        inner: Cursor<Vec<u8>>,
        seeks: u32,
    }

    impl Read for CountingCursor {
        fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
            self.inner.read(buf)
        }
    }

    impl Seek for CountingCursor {
        fn seek(&mut self, pos: SeekFrom) -> std::io::Result<u64> {
            self.seeks += 1;
            self.inner.seek(pos)
        }
    }

    fn counting(len: usize) -> CountingCursor {
        CountingCursor {
            inner: Cursor::new(vec![7u8; len]),
            seeks: 0,
        }
    }

    #[test]
    fn a_known_start_and_length_save_every_seek() {
        // What a per-scan reader pays for: each bounded read asks the stream
        // how long it is and where it is, which is a seek either way.
        let mut asked = BinaryReader::new(counting(64));
        assert_eq!(asked.read_bytes(4).unwrap().len(), 4);
        assert!(asked.inner.seeks > 0, "new() has to ask the stream");
        assert_eq!(asked.length().unwrap(), 64);

        let mut told = BinaryReader::at(counting(64), 0, 64);
        assert_eq!(told.length().unwrap(), 64);
        assert_eq!(told.read_bytes(4).unwrap().len(), 4);
        assert_eq!(told.position(), 4);
        assert_eq!(told.inner.seeks, 0, "a known start and length need none");
        // The cap on untrusted lengths still applies, and still counts from
        // where this reader started.
        assert!(told.read_bytes(1024).is_err());
        assert_eq!(told.remaining().unwrap(), 60);
    }

    #[test]
    fn a_reader_part_way_into_the_stream_knows_what_is_left() {
        let mut source = counting(64);
        source.seek(SeekFrom::Start(48)).unwrap();
        let mut r = BinaryReader::at(source, 48, 64);
        assert_eq!(r.remaining().unwrap(), 16);
        assert!(r.read_bytes(24).is_err(), "only 16 bytes are left");
        assert_eq!(r.read_bytes(16).unwrap().len(), 16);
        assert_eq!(r.remaining().unwrap(), 0);
        // Seeking makes the recorded start meaningless, so it is dropped and
        // the stream is asked again.
        r.seek_to(32).unwrap();
        assert_eq!(r.remaining().unwrap(), 32);
    }

    #[test]
    fn items_are_read_in_one_go_and_the_buffer_is_not_cleared() {
        let mut r = BinaryReader::at(counting(64), 0, 64);
        let mut buf = Vec::new();
        assert_eq!(r.read_items_into(&mut buf, 4, 8).unwrap(), &[7u8; 32]);
        assert_eq!(r.position(), 32);
        let capacity = buf.capacity();
        // A smaller run writes over the front and returns only its own
        // bytes; the buffer keeps its length, so nothing is zeroed again.
        assert_eq!(r.read_items_into(&mut buf, 2, 8).unwrap().len(), 16);
        assert_eq!(buf.len(), 32);
        assert_eq!(buf.capacity(), capacity, "the buffer should be reused");
        assert_eq!(r.inner.seeks, 0);
    }

    #[test]
    fn a_stream_that_ends_early_reports_the_item_it_could_not_finish() {
        // Three whole 8-byte items and half of a fourth. A loop reading one
        // item at a time would fail at offset 24 needing 8 bytes; so does
        // this.
        let mut r = BinaryReader::at(counting(28), 0, 32);
        let mut buf = Vec::new();
        let error = r.read_items_into(&mut buf, 4, 8).unwrap_err();
        assert!(
            matches!(
                error,
                Error::UnexpectedEof {
                    offset: 24,
                    needed: 8
                }
            ),
            "got {error:?}"
        );
        assert_eq!(r.position(), 24);
    }

    #[test]
    fn an_item_count_beyond_the_stream_is_capped_before_allocating() {
        let mut r = BinaryReader::at(counting(64), 0, 64);
        let mut buf = Vec::new();
        let error = r.read_items_into(&mut buf, 1_000_000, 8).unwrap_err();
        assert!(
            matches!(error, Error::AllocationTooLarge { .. }),
            "{error:?}"
        );
        assert!(buf.is_empty());
    }

    #[test]
    fn windows_filetime_zero_is_zero() {
        let bytes = 0u64.to_le_bytes().to_vec();
        let mut r = BinaryReader::new(Cursor::new(bytes));
        assert_eq!(r.read_windows_filetime().unwrap(), 0.0);
    }
}
