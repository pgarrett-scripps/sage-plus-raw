//! Strict decoding for the qualified v66 Orbitrap layout used by MSConvert RS.
use crate::error::{Error, Result};
use crate::reader::BinaryReader;
use crate::scan_event::{FractionCollector, Reaction, ScanEvent, ScanEventPreamble};
use std::io::{Read, Seek};

fn invalid(message: &str) -> Error {
    Error::InvalidScanEvent(message.into())
}

#[cfg(test)]
mod multiplex_tests {
    use super::*;
    #[test]
    fn preserves_each_reaction_and_rejects_truncation() {
        let mut p = vec![0u8; 136];
        p[6] = 2;
        let mut b = vec![0u8; 312];
        b[..4].copy_from_slice(&4u32.to_le_bytes());
        for (i, mz) in [522.987607f64, 698.067188, 502.978512, 768.099021]
            .into_iter()
            .enumerate()
        {
            let offset = 4 + 56 * i;
            b[offset..offset + 8].copy_from_slice(&mz.to_le_bytes());
            b[offset + 8..offset + 16].copy_from_slice(&5f64.to_le_bytes());
            b[offset + 16..offset + 24].copy_from_slice(&27f64.to_le_bytes());
        }
        b[228..232].copy_from_slice(&1u32.to_le_bytes());
        b[232..240].copy_from_slice(&100f64.to_le_bytes());
        b[240..248].copy_from_slice(&1000f64.to_le_bytes());
        b[248..252].copy_from_slice(&5u32.to_le_bytes());
        b[292..296].copy_from_slice(&1u32.to_le_bytes());
        let event = multiplexed(ScanEventPreamble { bytes: p.clone() }, &b, 4).unwrap();
        assert_eq!(event.reactions.len(), 4);
        assert_eq!(event.reactions[3].precursor_mz, 768.099021);
        assert!(multiplexed(ScanEventPreamble { bytes: p }, &b[..311], 4).is_err());
    }
}

/// Qualified v66 MSX records have 56-byte reactions followed by one mass range.
pub(crate) fn multiplexed(
    preamble: ScanEventPreamble,
    body: &[u8],
    count: usize,
) -> Result<ScanEvent> {
    let fc = count
        .checked_mul(56)
        .and_then(|n| n.checked_add(8))
        .ok_or_else(|| invalid("reaction count overflow"))?;
    if preamble.bytes.get(6).copied() != Some(2) || u32_at(body, fc - 4)? != 1 {
        return Err(invalid("unsupported multiplexed event layout"));
    }
    let low_mz = f64_at(body, fc)?;
    let high_mz = f64_at(body, fc + 8)?;
    let n = u32_at(body, fc + 16)? as usize;
    if low_mz <= 0.0
        || high_mz < low_mz
        || high_mz > 50_000.0
        || n != 5
        || body.len() != fc + 80
        || u32_at(body, fc + 60)? != 1
        || body[fc + 64..].iter().any(|&b| b != 0)
    {
        return Err(invalid("unqualified multiplexed range or calibration"));
    }
    let mut reactions = Vec::with_capacity(count);
    for i in 0..count {
        let offset = 4 + 56 * i;
        let precursor_mz = f64_at(body, offset)?;
        let width = f64_at(body, offset + 8)?;
        let energy = f64_at(body, offset + 16)?;
        if precursor_mz <= 0.0
            || precursor_mz > 50_000.0
            || width <= 0.0
            || energy < 0.0
            || body[offset + 32..offset + 56].iter().any(|&b| b != 0)
        {
            return Err(invalid("invalid multiplexed reaction"));
        }
        reactions.push(Reaction {
            precursor_mz,
            unknown_double: width,
            energy,
            unknown_long1: u32_at(body, offset + 24)?,
            unknown_long2: u32_at(body, offset + 28)?,
        });
    }
    let coefficients = (0..n)
        .map(|i| f64_at(body, fc + 20 + 8 * i))
        .collect::<Result<Vec<_>>>()?;
    Ok(ScanEvent {
        preamble,
        reactions,
        fraction_collectors: vec![FractionCollector { low_mz, high_mz }],
        coefficients,
    })
}

fn u32_at(bytes: &[u8], offset: usize) -> Result<u32> {
    let bytes = bytes
        .get(offset..offset + 4)
        .ok_or_else(|| invalid("truncated integer"))?;
    Ok(u32::from_le_bytes(
        bytes.try_into().map_err(|_| invalid("truncated integer"))?,
    ))
}

fn f64_at(bytes: &[u8], offset: usize) -> Result<f64> {
    let bytes = bytes
        .get(offset..offset + 8)
        .ok_or_else(|| invalid("truncated float"))?;
    let value = f64::from_le_bytes(bytes.try_into().map_err(|_| invalid("truncated float"))?);
    if !value.is_finite() {
        return Err(invalid("nonfinite event field"));
    }
    Ok(value)
}

/// Probe the complete record signature before selecting this layout.
/// Once selected, every record and the stream boundary must validate.
pub(crate) fn read<R: Read + Seek>(
    reader: &mut BinaryReader<R>,
    version: u32,
    count: u32,
    end: u64,
) -> Result<Option<Vec<ScanEvent>>> {
    if version != 66 || count == 0 {
        return Ok(None);
    }
    let start = reader.position();
    let length = end
        .checked_sub(start)
        .ok_or_else(|| invalid("reversed event bounds"))?;
    if length < u64::from(count) * 232 || length > u64::from(count) * 288 {
        return Ok(None);
    }
    let first = reader.read_bytes(164)?;
    reader.seek_to(start)?;
    if first[6] != 1
        || first[10] != 0
        || first[40] != 4
        || u32_at(&first, 136)? != 0
        || u32_at(&first, 140)? != 1
        || u32_at(&first, 160)? != 7
    {
        return Ok(None);
    }
    let mut events = Vec::new();
    for _ in 0..count {
        if end.saturating_sub(reader.position()) < 136 {
            return Err(invalid("truncated preamble"));
        }
        let preamble = reader.read_bytes(136)?;
        if !matches!(preamble[6], 1 | 2)
            || preamble[4] > 1
            || preamble[5] > 1
            || preamble[40] != 4
            || preamble[10] != u8::from(preamble[6] == 2)
        {
            return Err(invalid("unsupported qualified event preamble"));
        }
        let size = if preamble[6] == 2 { 152 } else { 96 };
        if end.saturating_sub(reader.position()) < size as u64 {
            return Err(invalid("truncated event body"));
        }
        let body = reader.read_bytes(size)?;
        events.push(decode(preamble, &body)?);
    }
    if reader.position() != end {
        return Err(invalid("unconsumed event bytes"));
    }
    Ok(Some(events))
}

fn decode(preamble: Vec<u8>, body: &[u8]) -> Result<ScanEvent> {
    let ms2 = preamble[6] == 2;
    let fc = if ms2 { 64 } else { 8 };
    if body.len() != if ms2 { 152 } else { 96 }
        || u32_at(body, 0)? != u32::from(ms2)
        || u32_at(body, fc - 4)? != 1
        || u32_at(body, fc + 16)? != 7
    {
        return Err(invalid("unsupported record fields"));
    }
    let low_mz = f64_at(body, fc)?;
    let high_mz = f64_at(body, fc + 8)?;
    if low_mz <= 0.0 || high_mz < low_mz || high_mz > 50_000.0 {
        return Err(invalid("invalid scan range"));
    }
    let coefficients = (0..7)
        .map(|i| f64_at(body, fc + 20 + i * 8))
        .collect::<Result<Vec<_>>>()?;
    if coefficients[3] <= 0.0 || body[fc + 76..].iter().any(|&b| b != 0) {
        return Err(invalid("invalid calibration or event tail"));
    }
    let mut reactions = Vec::new();
    if ms2 {
        let precursor_mz = f64_at(body, 4)?;
        let width = f64_at(body, 12)?;
        let energy = f64_at(body, 20)?;
        if !(0.0..=50_000.0).contains(&precursor_mz)
            || precursor_mz == 0.0
            || width < 0.0
            || energy < 0.0
            || body[36..60].iter().any(|&b| b != 0)
        {
            return Err(invalid("unsupported reaction metadata"));
        }
        reactions.push(Reaction {
            precursor_mz,
            unknown_double: width,
            energy,
            unknown_long1: u32_at(body, 28)?,
            unknown_long2: u32_at(body, 32)?,
        });
    }
    Ok(ScanEvent {
        preamble: ScanEventPreamble { bytes: preamble },
        reactions,
        fraction_collectors: vec![FractionCollector { low_mz, high_mz }],
        coefficients,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn event(ms2: bool) -> Vec<u8> {
        let mut bytes = vec![0u8; if ms2 { 288 } else { 232 }];
        bytes[6] = if ms2 { 2 } else { 1 };
        bytes[10] = u8::from(ms2);
        bytes[40] = 4;
        let body = &mut bytes[136..];
        body[..4].copy_from_slice(&u32::from(ms2).to_le_bytes());
        let fc = if ms2 { 64 } else { 8 };
        body[fc - 4..fc].copy_from_slice(&1u32.to_le_bytes());
        body[fc..fc + 8].copy_from_slice(&100f64.to_le_bytes());
        body[fc + 8..fc + 16].copy_from_slice(&1000f64.to_le_bytes());
        body[fc + 16..fc + 20].copy_from_slice(&7u32.to_le_bytes());
        body[fc + 44..fc + 52].copy_from_slice(&2e8f64.to_le_bytes());
        if ms2 {
            body[4..12].copy_from_slice(&500f64.to_le_bytes());
            body[12..20].copy_from_slice(&2f64.to_le_bytes());
            body[20..28].copy_from_slice(&30f64.to_le_bytes());
        }
        bytes
    }

    #[test]
    fn mixed_records_do_not_use_the_average_length() {
        let bytes = [event(false), event(true), event(false)].concat();
        let end = bytes.len() as u64;
        let mut r = BinaryReader::new(Cursor::new(bytes));
        let events = read(&mut r, 66, 3, end).unwrap().unwrap();
        assert_eq!(events[1].reactions[0].precursor_mz, 500.0);
        assert_eq!(events[1].reactions[0].energy, 30.0);
        assert_eq!(events[2].coefficients[3], 2e8);
        assert_eq!(r.position(), end);
    }

    #[test]
    fn malformed_selected_layout_does_not_fall_back() {
        let mut bytes = [event(false), event(true)].concat();
        bytes[232 + 136 + 36] = 1;
        let end = bytes.len() as u64;
        assert!(read(&mut BinaryReader::new(Cursor::new(bytes)), 66, 2, end).is_err());
    }

    #[test]
    fn rejects_truncation_and_trailing_data() {
        for change in [-1i64, 1] {
            let mut bytes = [event(false), event(true)].concat();
            bytes.resize((bytes.len() as i64 + change) as usize, 0);
            let end = bytes.len() as u64;
            assert!(read(&mut BinaryReader::new(Cursor::new(bytes)), 66, 2, end).is_err());
        }
    }
}
