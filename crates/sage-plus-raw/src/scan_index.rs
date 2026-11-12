use crate::binary_reader::BinaryReader;
use crate::error::Result;
use std::io::{Read, Seek};

/// One entry in the scan index array: 88 bytes in version 66, 80 in version
/// 64, 72 before that (`min_size_for_version`).
///
/// `RawFileReader::open` reads the entries back to back from
/// `RunHeader::scan_index_addr`, so what each version reads is also the
/// stride. Measured on real acquisitions: in every version 66 file the bytes
/// between `scan_index_addr` and the trailer stream that follows it are
/// exactly `num_scans * 88`, and in a version 64 file `num_scans * 80`.
#[derive(Debug)]
pub struct ScanIndexEntry {
    pub index: u32,
    pub scan_event: u16,
    pub scan_segment: u16,
    pub data_size: u32,
    pub start_time: f64,
    pub total_current: f64,
    pub base_intensity: f64,
    pub base_mz: f64,
    pub low_mz: f64,
    pub high_mz: f64,
    /// 64-bit offset relative to data stream start.
    pub offset: u64,
}

impl ScanIndexEntry {
    /// Minimum on-disk size in bytes of one entry for the given file version.
    /// Used to bound `Vec::with_capacity(num_scans)` against the file's
    /// actual remaining size before allocating.
    pub(crate) fn min_size_for_version(version: u32) -> u64 {
        // offset_32(4) + index(4) + scan_event(2) + scan_segment(2) + next(4)
        // + unk(4) + data_size(4) + 6 * f64(8) = 72 bytes, common to every
        // version (the offset_32 field is reused directly pre-v64).
        const BASE: u64 = 72;
        if version >= 66 {
            BASE + 8 /* u64 offset instead of offset_32 */ + 8 /* 2 trailing u32s */
        } else if version >= 64 {
            BASE + 8 /* u64 offset instead of offset_32 */
        } else {
            BASE
        }
    }

    pub(crate) fn read<R: Read + Seek>(r: &mut BinaryReader<R>, version: u32) -> Result<Self> {
        let offset_32 = r.read_u32()?;
        let index = r.read_u32()?;
        let scan_event = r.read_u16()?;
        let scan_segment = r.read_u16()?;
        let _next = r.read_u32()?;
        let _unk = r.read_u32()?;
        let data_size = r.read_u32()?;
        let start_time = r.read_f64()?;
        let total_current = r.read_f64()?;
        let base_intensity = r.read_f64()?;
        let base_mz = r.read_f64()?;
        let low_mz = r.read_f64()?;
        let high_mz = r.read_f64()?;

        let offset = if version >= 64 {
            r.read_u64()?
        } else {
            offset_32 as u64
        };

        if version >= 66 {
            // Two trailing unknown u32s
            let _unk1 = r.read_u32()?;
            let _unk2 = r.read_u32()?;
        }

        Ok(Self {
            index,
            scan_event,
            scan_segment,
            data_size,
            start_time,
            total_current,
            base_intensity,
            base_mz,
            low_mz,
            high_mz,
            offset,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    /// The bytes of one scan-index entry as `version` writes them.
    ///
    /// Common to every version: a 32-bit offset, the index, the event and
    /// segment words, two unused words, the data size, then six f64 fields.
    /// From v64 a 64-bit offset follows and supersedes the 32-bit one; v66
    /// adds two trailing words.
    fn entry_bytes(version: u32, offset_32: u32, offset_64: u64) -> Vec<u8> {
        let mut bytes = offset_32.to_le_bytes().to_vec();
        bytes.extend_from_slice(&7u32.to_le_bytes()); // index
        bytes.extend_from_slice(&2u16.to_le_bytes()); // scan_event
        bytes.extend_from_slice(&1u16.to_le_bytes()); // scan_segment
        bytes.extend_from_slice(&0u32.to_le_bytes()); // next
        bytes.extend_from_slice(&0u32.to_le_bytes()); // unk
        bytes.extend_from_slice(&1234u32.to_le_bytes()); // data_size
        for value in [12.5f64, 1.0e6, 5.0e5, 445.12, 100.0, 1500.0] {
            bytes.extend_from_slice(&value.to_le_bytes());
        }
        if version >= 64 {
            bytes.extend_from_slice(&offset_64.to_le_bytes());
        }
        if version >= 66 {
            bytes.extend_from_slice(&0u32.to_le_bytes());
            bytes.extend_from_slice(&0u32.to_le_bytes());
        }
        bytes
    }

    /// Every field but the offset is at the same place in every version.
    fn assert_common_fields(entry: &ScanIndexEntry) {
        assert_eq!(entry.index, 7);
        assert_eq!(entry.scan_event, 2);
        assert_eq!(entry.scan_segment, 1);
        assert_eq!(entry.data_size, 1234);
        assert_eq!(entry.start_time, 12.5);
        assert_eq!(entry.total_current, 1.0e6);
        assert_eq!(entry.base_intensity, 5.0e5);
        assert_eq!(entry.base_mz, 445.12);
        assert_eq!(entry.low_mz, 100.0);
        assert_eq!(entry.high_mz, 1500.0);
    }

    #[test]
    fn min_size_for_version_is_the_entry_stride_of_each_layout() {
        // These are the strides the files themselves show: the bytes between
        // `scan_index_addr` and the following trailer stream divided by the
        // scan count. 88 for version 66 (the six public ProteoWizard
        // fixtures, and a 67,599-scan Q Exactive acquisition: 5,948,712
        // bytes), 80 for version 64 (a 38,687-scan LTQ Orbitrap Velos
        // acquisition: 3,094,960 bytes).
        for version in [8, 47, 57, 60, 62, 63] {
            assert_eq!(
                ScanIndexEntry::min_size_for_version(version),
                72,
                "v{version}"
            );
        }
        assert_eq!(ScanIndexEntry::min_size_for_version(64), 80);
        assert_eq!(ScanIndexEntry::min_size_for_version(66), 88);
    }

    #[test]
    fn read_pre_v64_takes_the_offset_from_the_32_bit_field() {
        let bytes = entry_bytes(63, 4096, 0);
        let written = bytes.len() as u64;
        let mut r = BinaryReader::new(Cursor::new(bytes));
        let entry = ScanIndexEntry::read(&mut r, 63).unwrap();
        assert_common_fields(&entry);
        assert_eq!(entry.offset, 4096);
        assert_eq!(written, 72);
        // What the allocation bound in `RawFileReader::open` relies on: an
        // entry never takes fewer bytes than the bound assumes.
        assert!(r.position() >= ScanIndexEntry::min_size_for_version(63));
    }

    #[test]
    fn read_v64_takes_the_offset_from_the_64_bit_field() {
        // The 32-bit field is left at a value the 64-bit one must override,
        // so a reader that ignored the wider field would be caught.
        let bytes = entry_bytes(64, 4096, 8_589_934_592);
        let written = bytes.len() as u64;
        let mut r = BinaryReader::new(Cursor::new(bytes));
        let entry = ScanIndexEntry::read(&mut r, 64).unwrap();
        assert_common_fields(&entry);
        assert_eq!(entry.offset, 8_589_934_592);
        assert_eq!(written, 80);
        assert!(r.position() >= ScanIndexEntry::min_size_for_version(64));
    }

    #[test]
    fn read_v66_takes_the_64_bit_offset_and_consumes_the_trailing_words() {
        let bytes = entry_bytes(66, 4096, 8_589_934_592);
        let written = bytes.len() as u64;
        let mut r = BinaryReader::new(Cursor::new(bytes));
        let entry = ScanIndexEntry::read(&mut r, 66).unwrap();
        assert_common_fields(&entry);
        assert_eq!(entry.offset, 8_589_934_592);
        assert_eq!(written, 88);
        // The two trailing words are part of the entry: the next entry
        // starts after them, so the reader must have consumed them.
        assert_eq!(r.position(), written);
        assert!(r.position() >= ScanIndexEntry::min_size_for_version(66));
    }

    #[test]
    fn read_rejects_an_entry_that_ends_early() {
        let mut bytes = entry_bytes(66, 4096, 8_589_934_592);
        bytes.truncate(bytes.len() - 4); // one trailing word short
        let mut r = BinaryReader::new(Cursor::new(bytes));
        let error = ScanIndexEntry::read(&mut r, 66).unwrap_err();
        assert!(
            matches!(error, crate::Error::UnexpectedEof { .. }),
            "got {error:?}"
        );
    }
}
