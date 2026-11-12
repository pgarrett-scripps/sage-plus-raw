use crate::binary_reader::BinaryReader;
use crate::error::{Error, Result};
use crate::types::GenericType;
use std::io::{Cursor, Read, Seek};

/// A field descriptor within a GenericDataHeader.
#[derive(Debug, Clone)]
pub struct GenericDataDescriptor {
    pub field_type: GenericType,
    pub length: u32,
    pub label: String,
}

/// Self-describing header for GenericRecord streams.
#[derive(Debug)]
pub struct GenericDataHeader {
    pub fields: Vec<GenericDataDescriptor>,
}

/// A typed value from a generic record.
#[derive(Debug, Clone, PartialEq)]
pub enum GenericValue {
    Gap,
    Int8(i8),
    Bool(bool),
    UInt8(u8),
    Int16(i16),
    UInt16(u16),
    Int32(i32),
    UInt32(u32),
    Float32(f32),
    Float64(f64),
    String(String),
}

/// A single record decoded using a GenericDataHeader.
#[derive(Debug)]
pub struct GenericRecord {
    pub values: Vec<(String, GenericValue)>,
}

impl GenericDataHeader {
    /// Try to read a GenericDataHeader. Returns `None` if the data at the
    /// current position does not look like a valid header (e.g. unreasonable
    /// field count or invalid type codes). The reader position is restored
    /// on failure.
    pub(crate) fn try_read<R: Read + Seek>(r: &mut BinaryReader<R>) -> Result<Option<Self>> {
        let saved_pos = r.position();
        let n = r.read_u32()?;
        // A genuine schema has at least a couple of fields and fewer than ~500.
        // The error-log "gap" region that precedes the schema in v64+ can
        // otherwise mislead us by looking like a 0- or 1-field header.
        if !(2..=500).contains(&n) {
            r.seek_to(saved_pos)?;
            return Ok(None);
        }
        let mut fields = Vec::with_capacity(n as usize);
        for _ in 0..n {
            let type_code = r.read_u32()?;
            match GenericType::from_u32(type_code) {
                Some(field_type) => {
                    let length = r.read_u32()?;
                    // Character count of the label. Real Thermo labels are
                    // short and printable; require it to look sane or the
                    // whole header is bogus.
                    let label_start = r.position();
                    let char_count = r.read_u32()?;
                    if char_count > 200 {
                        r.seek_to(saved_pos)?;
                        return Ok(None);
                    }
                    r.seek_to(label_start)?;
                    let label = match r.read_pascal_string() {
                        Ok(s) => s,
                        Err(crate::error::Error::InvalidUtf16(_)) => {
                            r.seek_to(saved_pos)?;
                            return Ok(None);
                        }
                        Err(e) => return Err(e),
                    };
                    if !label_is_plausible(&label) {
                        r.seek_to(saved_pos)?;
                        return Ok(None);
                    }
                    fields.push(GenericDataDescriptor {
                        field_type,
                        length,
                        label,
                    });
                }
                None => {
                    r.seek_to(saved_pos)?;
                    return Ok(None);
                }
            }
        }
        let hdr = Self { fields };
        if !hdr.looks_meaningful() {
            r.seek_to(saved_pos)?;
            return Ok(None);
        }
        Ok(Some(hdr))
    }

    /// A schema is "meaningful" if it contains at least a few fields with
    /// real labels and has a non-trivial fixed record size. Used to reject
    /// false positives picked up by the forward scan.
    fn looks_meaningful(&self) -> bool {
        let named = self.fields.iter().filter(|f| !f.label.is_empty()).count();
        named >= 2 && self.fixed_record_size() > 0
    }

    /// Sum of fixed byte sizes contributed by each descriptor. For variable
    /// types (String/WideString) the descriptor's `length` field is used as
    /// the storage allocation - which is the fixed on-disk size per record.
    pub(crate) fn fixed_record_size(&self) -> usize {
        self.fields
            .iter()
            .map(|f| match f.field_type {
                GenericType::Gap => 0,
                GenericType::Int8
                | GenericType::Bool
                | GenericType::BoolYesNo
                | GenericType::BoolOnOff
                | GenericType::UInt8 => 1,
                GenericType::Int16 | GenericType::UInt16 => 2,
                GenericType::Int32 | GenericType::UInt32 | GenericType::Float32 => 4,
                GenericType::Float64 => 8,
                GenericType::AsciiString => f.length as usize,
                GenericType::WideString => f.length as usize * 2,
            })
            .sum()
    }

    /// Scan forward from the current position for a plausible GenericDataHeader
    /// in a bounded window. The v64+ error-log region contains padding bytes
    /// before the scan-parameters schema whose size isn't easily computed, so
    /// we locate the schema by scanning for a valid signature.
    pub(crate) fn find_forward<R: Read + Seek>(
        r: &mut BinaryReader<R>,
        max_scan: u64,
        expected_record_size: Option<usize>,
    ) -> Result<Option<Self>> {
        let start = r.position();
        let cap = max_scan.min(4 * 1024 * 1024) as usize;
        r.seek_to(start)?;
        let buf = r.read_bytes(cap)?;
        // Two passes: first require the schema's fixed record size to match
        // the tail; second accept any meaningful schema.
        //
        // Parse candidates entirely from the in-memory buffer using a Cursor so
        // that we never seek the underlying file reader for each false positive.
        // This avoids O(n) file seeks when the error-log gap is large (>1 MB).
        for pass in 0..2 {
            let mut offset = 0usize;
            while offset + 4 <= buf.len() {
                let n = crate::bytes::read_u32_le(&buf, offset)?;
                if (2..=500).contains(&n) {
                    let mut cursor = BinaryReader::new(Cursor::new(&buf[offset..]));
                    // A candidate that runs past the end of the window is not
                    // a header: skip it rather than fail the whole search.
                    let candidate = match Self::try_read(&mut cursor) {
                        Err(Error::UnexpectedEof { .. } | Error::AllocationTooLarge { .. }) => None,
                        other => other?,
                    };
                    if let Some(hdr) = candidate {
                        let size_ok = match (pass, expected_record_size) {
                            (0, Some(want)) => hdr.fixed_record_size() == want,
                            _ => true,
                        };
                        if size_ok {
                            return Ok(Some(hdr));
                        }
                    }
                }
                offset += 2;
            }
            if expected_record_size.is_none() {
                break;
            }
        }
        r.seek_to(start)?;
        Ok(None)
    }
}

/// Heuristic: a GDH field label must either be empty or have reasonable
/// length. Labels are sometimes short single-character sentinels so we
/// don't require printability.
fn label_is_plausible(s: &str) -> bool {
    s.len() <= 200
}

/// Decode one field over the top of `slot`, reusing the `String` already
/// there for a string field and `scratch` for its bytes. The values and the
/// errors are those [`GenericRecord::read`] has always produced; this is the
/// one place a generic-record field is decoded.
fn read_field_into<R: Read + Seek>(
    r: &mut BinaryReader<R>,
    desc: &GenericDataDescriptor,
    slot: &mut GenericValue,
    scratch: &mut Vec<u8>,
) -> Result<()> {
    let value = match desc.field_type {
        GenericType::Gap => GenericValue::Gap,
        GenericType::Int8 => GenericValue::Int8(r.read_i8()?),
        GenericType::Bool | GenericType::BoolYesNo | GenericType::BoolOnOff => {
            GenericValue::Bool(r.read_u8()? != 0)
        }
        GenericType::UInt8 => GenericValue::UInt8(r.read_u8()?),
        GenericType::Int16 => GenericValue::Int16(r.read_i16()?),
        GenericType::UInt16 => GenericValue::UInt16(r.read_u16()?),
        GenericType::Int32 => GenericValue::Int32(r.read_i32()?),
        GenericType::UInt32 => GenericValue::UInt32(r.read_u32()?),
        GenericType::Float32 => GenericValue::Float32(r.read_f32()?),
        GenericType::Float64 => GenericValue::Float64(r.read_f64()?),
        GenericType::AsciiString => {
            if desc.length > 0 {
                let bytes = r.read_bytes_scratch(scratch, desc.length as usize)?;
                let end = bytes.iter().position(|&b| b == 0).unwrap_or(bytes.len());
                let text = String::from_utf8_lossy(&bytes[..end]);
                store_text(slot, &text);
            } else {
                store_text(slot, "");
            }
            return Ok(());
        }
        GenericType::WideString => {
            if desc.length > 0 {
                // Lossy: a malformed trailer string must not fail the whole
                // open.
                let mut text = take_text(slot);
                let result =
                    r.read_utf16_fixed_lossy_into(desc.length as usize * 2, scratch, &mut text);
                *slot = GenericValue::String(text);
                result?;
            } else {
                store_text(slot, "");
            }
            return Ok(());
        }
    };
    *slot = value;
    Ok(())
}

/// Put `text` in `slot`, keeping the `String` already there if it holds one.
fn store_text(slot: &mut GenericValue, text: &str) {
    match slot {
        GenericValue::String(existing) => {
            existing.clear();
            existing.push_str(text);
        }
        other => *other = GenericValue::String(text.to_string()),
    }
}

/// The `String` in `slot`, to be written and put back, or a new one.
fn take_text(slot: &mut GenericValue) -> String {
    match std::mem::replace(slot, GenericValue::Gap) {
        GenericValue::String(existing) => existing,
        _ => String::new(),
    }
}

impl GenericRecord {
    pub(crate) fn read<R: Read + Seek>(
        r: &mut BinaryReader<R>,
        header: &GenericDataHeader,
    ) -> Result<Self> {
        let mut record = Self { values: Vec::new() };
        record.read_into(r, header, &mut Vec::new())?;
        Ok(record)
    }

    /// Read one record over the top of this one, for a caller walking a
    /// whole stream: a per-scan stream of 136k records then allocates
    /// nothing per record instead of this record's `Vec`, one `String` per
    /// field for the label, and one per string value.
    ///
    /// The labels are cloned from `header` only when this record does not
    /// already carry exactly them, so passing the same header every time
    /// costs one comparison per field. `scratch` is a byte buffer for
    /// fixed-length string fields; its contents are not meaningful
    /// afterwards. The values are those [`Self::read`] produces.
    pub(crate) fn read_into<R: Read + Seek>(
        &mut self,
        r: &mut BinaryReader<R>,
        header: &GenericDataHeader,
        scratch: &mut Vec<u8>,
    ) -> Result<()> {
        if !self.describes(header) {
            self.values.clear();
            self.values.reserve(header.fields.len());
            for desc in &header.fields {
                self.values.push((desc.label.clone(), GenericValue::Gap));
            }
        }
        for (slot, desc) in self.values.iter_mut().zip(&header.fields) {
            read_field_into(r, desc, &mut slot.1, scratch)?;
        }
        Ok(())
    }

    /// Whether this record's labels are already `header`'s, field for field.
    fn describes(&self, header: &GenericDataHeader) -> bool {
        self.values.len() == header.fields.len()
            && self
                .values
                .iter()
                .zip(&header.fields)
                .all(|((label, _), desc)| *label == desc.label)
    }

    /// Look up a field by label and return a reference to its value.
    pub fn get(&self, label: &str) -> Option<&GenericValue> {
        self.values.iter().find(|(l, _)| l == label).map(|(_, v)| v)
    }

    /// Get a float64 field by label.
    pub fn get_f64(&self, label: &str) -> Option<f64> {
        match self.get(label)? {
            GenericValue::Float64(v) => Some(*v),
            GenericValue::Float32(v) => Some(*v as f64),
            _ => None,
        }
    }

    /// Get a float32 field by label.
    pub fn get_f32(&self, label: &str) -> Option<f32> {
        match self.get(label)? {
            GenericValue::Float32(v) => Some(*v),
            GenericValue::Float64(v) => Some(*v as f32),
            _ => None,
        }
    }

    /// Get an i32 field by label.
    pub fn get_i32(&self, label: &str) -> Option<i32> {
        match self.get(label)? {
            GenericValue::Int32(v) => Some(*v),
            GenericValue::Int16(v) => Some(*v as i32),
            GenericValue::Int8(v) => Some(*v as i32),
            _ => None,
        }
    }

    /// Get a string field by label.
    pub fn get_string(&self, label: &str) -> Option<&str> {
        match self.get(label)? {
            GenericValue::String(v) => Some(v.as_str()),
            _ => None,
        }
    }
}

impl GenericValue {
    /// Get as f64, converting numeric types.
    pub fn as_f64(&self) -> Option<f64> {
        match self {
            Self::Float64(v) => Some(*v),
            Self::Float32(v) => Some(*v as f64),
            Self::Int32(v) => Some(*v as f64),
            Self::UInt32(v) => Some(*v as f64),
            Self::Int16(v) => Some(*v as f64),
            Self::UInt16(v) => Some(*v as f64),
            Self::Int8(v) => Some(*v as f64),
            Self::UInt8(v) => Some(*v as f64),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_util::pascal_string;
    use std::io::Cursor;

    /// One field descriptor: type_code(u32) + length(u32) + label (pascal string).
    fn field_bytes(type_code: u32, length: u32, label: &str) -> Vec<u8> {
        let mut out = type_code.to_le_bytes().to_vec();
        out.extend_from_slice(&length.to_le_bytes());
        out.extend_from_slice(&pascal_string(label));
        out
    }

    /// A full GDH blob: field count(u32) + that many field descriptors.
    fn gdh_bytes(fields: &[(u32, u32, &str)]) -> Vec<u8> {
        let mut out = (fields.len() as u32).to_le_bytes().to_vec();
        for &(t, l, label) in fields {
            out.extend_from_slice(&field_bytes(t, l, label));
        }
        out
    }

    #[test]
    fn wide_string_with_unpaired_surrogate_decodes_lossily() {
        let header = GenericDataHeader {
            fields: vec![
                GenericDataDescriptor {
                    field_type: GenericType::WideString,
                    length: 4,
                    label: "Comment:".into(),
                },
                GenericDataDescriptor {
                    field_type: GenericType::Int32,
                    length: 0,
                    label: "Next:".into(),
                },
            ],
        };
        // "A", lone high surrogate, "B", NUL terminator; then an i32.
        let mut bytes = Vec::new();
        for u in [0x0041u16, 0xD800, 0x0042, 0x0000] {
            bytes.extend_from_slice(&u.to_le_bytes());
        }
        bytes.extend_from_slice(&7i32.to_le_bytes());
        let mut r = BinaryReader::new(Cursor::new(bytes));
        let rec = GenericRecord::read(&mut r, &header).expect("lossy decode");
        match rec.get("Comment:") {
            Some(GenericValue::String(s)) => assert_eq!(s, "A\u{FFFD}B"),
            other => panic!("unexpected value: {other:?}"),
        }
        // The stream stays aligned for the following field.
        assert!(matches!(rec.get("Next:"), Some(GenericValue::Int32(7))));
    }

    /// A schema with one field of every kind a trailer stream mixes, so a
    /// record exercises both string paths and the scalars between them.
    fn mixed_header() -> GenericDataHeader {
        GenericDataHeader {
            fields: vec![
                GenericDataDescriptor {
                    field_type: GenericType::Float64,
                    length: 0,
                    label: "RT:".into(),
                },
                GenericDataDescriptor {
                    field_type: GenericType::AsciiString,
                    length: 6,
                    label: "Tag:".into(),
                },
                GenericDataDescriptor {
                    field_type: GenericType::WideString,
                    length: 4,
                    label: "Comment:".into(),
                },
                GenericDataDescriptor {
                    field_type: GenericType::Int32,
                    length: 0,
                    label: "Scan:".into(),
                },
            ],
        }
    }

    /// One record for `mixed_header`: an f64, a 6-byte ASCII field, four
    /// UTF-16 code units, an i32.
    fn mixed_record_bytes(rt: f64, tag: &str, comment: &str, scan: i32) -> Vec<u8> {
        let mut out = rt.to_le_bytes().to_vec();
        let mut ascii = tag.as_bytes().to_vec();
        ascii.resize(6, 0);
        out.extend_from_slice(&ascii);
        out.extend_from_slice(&crate::test_util::utf16_fixed(comment, 8));
        out.extend_from_slice(&scan.to_le_bytes());
        out
    }

    #[test]
    fn a_reused_record_decodes_what_a_fresh_one_does() {
        let header = mixed_header();
        let records = [
            mixed_record_bytes(1.5, "abc", "hi", 7),
            mixed_record_bytes(2.5, "de", "there", 8),
            // Shorter values than the record before, to catch a reused
            // `String` keeping the tail of the previous one.
            mixed_record_bytes(3.5, "f", "x", 9),
        ];
        let stream: Vec<u8> = records.concat();

        let mut fresh = BinaryReader::new(Cursor::new(stream.clone()));
        let one_by_one: Vec<GenericRecord> = (0..records.len())
            .map(|_| GenericRecord::read(&mut fresh, &header).expect("fresh record"))
            .collect();

        let mut reused = GenericRecord { values: Vec::new() };
        let mut scratch = Vec::new();
        let mut r = BinaryReader::new(Cursor::new(stream));
        let mut label_ptr = None;
        for expected in &one_by_one {
            reused
                .read_into(&mut r, &header, &mut scratch)
                .expect("reused record");
            assert_eq!(reused.values, expected.values);
            // Nothing is allocated per record after the first: the label
            // strings and the string value's buffer stay where they were.
            let ptr = (
                reused.values[0].0.as_ptr(),
                match &reused.values[1].1 {
                    GenericValue::String(s) => s.as_ptr(),
                    other => panic!("expected a string value: {other:?}"),
                },
            );
            if let Some(first) = label_ptr {
                assert_eq!(ptr, first);
            }
            label_ptr = Some(ptr);
        }
    }

    #[test]
    fn a_reused_record_rebuilds_itself_for_a_different_schema() {
        let header = mixed_header();
        let mut other = GenericDataHeader {
            fields: header.fields.clone(),
        };
        other.fields[0].label = "Start:".into();
        let bytes = mixed_record_bytes(1.5, "abc", "hi", 7);
        let mut reused = GenericRecord { values: Vec::new() };
        let mut scratch = Vec::new();
        let mut r = BinaryReader::new(Cursor::new(bytes.clone()));
        reused.read_into(&mut r, &header, &mut scratch).unwrap();
        assert!(reused.get("RT:").is_some());
        let mut r = BinaryReader::new(Cursor::new(bytes));
        reused.read_into(&mut r, &other, &mut scratch).unwrap();
        assert!(reused.get("RT:").is_none());
        assert!(reused.get("Start:").is_some());
    }

    fn try_read(bytes: Vec<u8>) -> Result<Option<GenericDataHeader>> {
        let mut r = BinaryReader::new(Cursor::new(bytes));
        GenericDataHeader::try_read(&mut r)
    }

    #[test]
    fn valid_two_field_header_parses() {
        let bytes = gdh_bytes(&[
            (GenericType::Float64 as u32, 0, "RT:"),
            (GenericType::Int32 as u32, 0, "Scan:"),
        ]);
        let hdr = try_read(bytes).unwrap().expect("should parse");
        assert_eq!(hdr.fields.len(), 2);
        assert_eq!(hdr.fields[0].label, "RT:");
        assert_eq!(hdr.fields[1].label, "Scan:");
        assert_eq!(hdr.fixed_record_size(), 8 + 4);
    }

    #[test]
    fn field_count_zero_is_rejected() {
        let bytes = gdh_bytes(&[]);
        assert!(try_read(bytes).unwrap().is_none());
    }

    #[test]
    fn field_count_one_is_rejected() {
        let bytes = gdh_bytes(&[(GenericType::Int32 as u32, 0, "Only:")]);
        assert!(try_read(bytes).unwrap().is_none());
    }

    #[test]
    fn field_count_above_500_is_rejected() {
        let mut bytes = 501u32.to_le_bytes().to_vec();
        bytes.extend_from_slice(&field_bytes(GenericType::Int32 as u32, 0, "A:"));
        assert!(try_read(bytes).unwrap().is_none());
    }

    #[test]
    fn invalid_type_code_is_rejected() {
        let bytes = gdh_bytes(&[(0xDEAD_BEEF, 0, "A:"), (GenericType::Int32 as u32, 0, "B:")]);
        assert!(try_read(bytes).unwrap().is_none());
    }

    #[test]
    fn label_char_count_above_200_is_rejected() {
        let mut out = 2u32.to_le_bytes().to_vec();
        out.extend_from_slice(&(GenericType::Int32 as u32).to_le_bytes());
        out.extend_from_slice(&0u32.to_le_bytes()); // length
        out.extend_from_slice(&300u32.to_le_bytes()); // implausible char count
        out.extend_from_slice(&[0u8; 20]); // some filler bytes (not enough for 300 chars)
        assert!(try_read(out).unwrap().is_none());
    }

    #[test]
    fn all_empty_labels_are_not_meaningful() {
        // named() requires at least 2 non-empty labels; both empty here.
        let bytes = gdh_bytes(&[
            (GenericType::Int32 as u32, 0, ""),
            (GenericType::Int32 as u32, 0, ""),
        ]);
        assert!(try_read(bytes).unwrap().is_none());
    }

    #[test]
    fn position_restored_on_rejection() {
        let mut bytes = vec![0xAAu8; 8]; // leading junk before the candidate
        let candidate_start = bytes.len();
        bytes.extend_from_slice(&gdh_bytes(&[])); // n=0, rejected
        let mut r = BinaryReader::new(Cursor::new(bytes));
        r.seek_to(candidate_start as u64).unwrap();
        let saved = r.position();
        assert!(GenericDataHeader::try_read(&mut r).unwrap().is_none());
        assert_eq!(r.position(), saved);
    }

    #[test]
    fn truncated_header_is_eof_not_panic() {
        let mut bytes = gdh_bytes(&[
            (GenericType::Float64 as u32, 0, "RT:"),
            (GenericType::Int32 as u32, 0, "Scan:"),
        ]);
        bytes.truncate(bytes.len() - 2);
        // Either a clean error or a `None` (truncated pascal string content
        // fails UTF-16 decoding / hits EOF) is acceptable; a panic is not.
        let _ = try_read(bytes);
    }

    #[test]
    fn fixed_record_size_sums_all_field_kinds() {
        let bytes = gdh_bytes(&[
            (GenericType::Int8 as u32, 0, "a:"),
            (GenericType::Int32 as u32, 0, "b:"),
            (GenericType::Float64 as u32, 0, "c:"),
            (GenericType::AsciiString as u32, 10, "d:"),
            (GenericType::WideString as u32, 5, "e:"),
        ]);
        let hdr = try_read(bytes).unwrap().unwrap();
        // 1 (Int8) + 4 (Int32) + 8 (Float64) + 10 (Ascii len) + 5*2 (Wide len*2)
        assert_eq!(hdr.fixed_record_size(), 1 + 4 + 8 + 10 + 10);
    }

    #[test]
    fn find_forward_locates_header_after_garbage_prefix() {
        let mut bytes = vec![0x11u8; 64]; // garbage that doesn't look like a valid GDH
        bytes.extend_from_slice(&gdh_bytes(&[
            (GenericType::Float64 as u32, 0, "RT:"),
            (GenericType::Int32 as u32, 0, "Scan:"),
        ]));
        let max_scan = bytes.len() as u64;
        let mut r = BinaryReader::new(Cursor::new(bytes));
        let hdr = GenericDataHeader::find_forward(&mut r, max_scan, None)
            .unwrap()
            .expect("should find the embedded header");
        assert_eq!(hdr.fields.len(), 2);
    }

    #[test]
    fn find_forward_skips_candidate_truncated_at_window_end() {
        // The real header does not match the expected record size, so the
        // first pass scans on to the end of the window...
        let mut bytes = gdh_bytes(&[
            (GenericType::Float64 as u32, 0, "RT:"),
            (GenericType::Int32 as u32, 0, "Scan:"),
        ]);
        // ...where a plausible field count is cut off mid-descriptor.
        bytes.extend_from_slice(&64u32.to_le_bytes());
        bytes.extend_from_slice(&(GenericType::Int32 as u32).to_le_bytes());
        bytes.extend_from_slice(&[0u8; 5]);
        let max_scan = bytes.len() as u64;
        let mut r = BinaryReader::new(Cursor::new(bytes));
        let hdr = GenericDataHeader::find_forward(&mut r, max_scan, Some(999))
            .unwrap()
            .expect("the second pass should still find the real header");
        assert_eq!(hdr.fields.len(), 2);
    }

    #[test]
    fn find_forward_returns_none_when_absent() {
        let bytes = vec![0x11u8; 128];
        let max_scan = bytes.len() as u64;
        let mut r = BinaryReader::new(Cursor::new(bytes));
        assert!(GenericDataHeader::find_forward(&mut r, max_scan, None)
            .unwrap()
            .is_none());
    }

    #[test]
    fn generic_record_reads_typed_fields() {
        let hdr = GenericDataHeader {
            fields: vec![
                GenericDataDescriptor {
                    field_type: GenericType::Int32,
                    length: 0,
                    label: "Scan:".to_string(),
                },
                GenericDataDescriptor {
                    field_type: GenericType::Float64,
                    length: 0,
                    label: "RT:".to_string(),
                },
                GenericDataDescriptor {
                    field_type: GenericType::AsciiString,
                    length: 8,
                    label: "Note:".to_string(),
                },
            ],
        };
        let mut bytes = 42i32.to_le_bytes().to_vec();
        bytes.extend_from_slice(&1.5f64.to_le_bytes());
        let mut note = b"hi".to_vec();
        note.resize(8, 0); // null-padded fixed-width ASCII field
        bytes.extend_from_slice(&note);

        let mut r = BinaryReader::new(Cursor::new(bytes));
        let record = GenericRecord::read(&mut r, &hdr).unwrap();
        assert_eq!(record.get_i32("Scan:"), Some(42));
        assert_eq!(record.get_f64("RT:"), Some(1.5));
        assert_eq!(record.get_string("Note:"), Some("hi"));
    }

    #[test]
    fn generic_record_truncated_is_eof() {
        let hdr = GenericDataHeader {
            fields: vec![GenericDataDescriptor {
                field_type: GenericType::Float64,
                length: 0,
                label: "RT:".to_string(),
            }],
        };
        let bytes = vec![0u8; 4]; // Float64 needs 8 bytes
        let mut r = BinaryReader::new(Cursor::new(bytes));
        assert!(GenericRecord::read(&mut r, &hdr).is_err());
    }

    #[test]
    fn generic_value_as_f64_converts_numeric_variants() {
        assert_eq!(GenericValue::Float64(2.5).as_f64(), Some(2.5));
        assert_eq!(GenericValue::Int32(7).as_f64(), Some(7.0));
        assert_eq!(GenericValue::UInt8(3).as_f64(), Some(3.0));
        assert_eq!(GenericValue::Gap.as_f64(), None);
        assert_eq!(GenericValue::String("x".into()).as_f64(), None);
    }
}
