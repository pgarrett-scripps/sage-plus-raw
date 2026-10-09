//! Self-describing decoding of the per-scan event stream.
//!
//! Every scan event on disk is a fixed preamble followed by counted lists:
//!
//! ```text
//! preamble                      (41 / 80 / 120 / 128 / 136 bytes by version)
//! u32 n_reactions,   n_reactions   x reaction     (32 bytes; 56 bytes in v66)
//! u32 n_mass_ranges, n_mass_ranges x (f64 low, f64 high)
//! u32 n_calibrators, n_calibrators x f64          (m/z conversion coefficients)
//! u32 n_source_frag, n_source_frag x f64          (source fragmentation values)
//! u32 n_sf_ranges,   n_sf_ranges   x (f64, f64)   (source fragmentation ranges)
//! v66 only: u32 n_chars, n_chars x u16           (UTF-16 text, e.g. a target name)
//! ```
//!
//! Event sizes therefore vary scan by scan (EThcD and SPS-MS3 events carry
//! two reactions, ion-trap events carry no calibrators, MSX events carry up
//! to dozens of reactions). The layout is accepted only when the whole stream
//! decodes and ends exactly at the next stream's address. Otherwise the
//! caller falls back to the legacy fixed-size heuristics.
use crate::scan_event::{FractionCollector, Reaction, ScanEvent, ScanEventPreamble};

/// Largest plausible list length in one event; guards against reading a
/// garbage count as a list that happens to fit the remaining bytes.
const MAX_LIST: usize = 4096;

struct Cursor<'a> {
    bytes: &'a [u8],
    pos: usize,
}

impl<'a> Cursor<'a> {
    fn take(&mut self, n: usize) -> Option<&'a [u8]> {
        let end = self.pos.checked_add(n)?;
        let slice = self.bytes.get(self.pos..end)?;
        self.pos = end;
        Some(slice)
    }

    fn u32(&mut self) -> Option<u32> {
        self.take(4)
            .map(|b| u32::from_le_bytes([b[0], b[1], b[2], b[3]]))
    }

    fn count(&mut self, item: usize) -> Option<usize> {
        let n = self.u32()? as usize;
        let remaining = self.bytes.len() - self.pos;
        (n <= MAX_LIST && n.checked_mul(item)? <= remaining).then_some(n)
    }
}

fn f64_at(bytes: &[u8], offset: usize) -> f64 {
    let mut b = [0u8; 8];
    b.copy_from_slice(&bytes[offset..offset + 8]);
    f64::from_le_bytes(b)
}

fn u32_at(bytes: &[u8], offset: usize) -> u32 {
    u32::from_le_bytes([
        bytes[offset],
        bytes[offset + 1],
        bytes[offset + 2],
        bytes[offset + 3],
    ])
}

fn decode_event(
    c: &mut Cursor<'_>,
    preamble_size: usize,
    reaction_size: usize,
    text: bool,
) -> Option<ScanEvent> {
    let preamble = c.take(preamble_size)?.to_vec();
    let n = c.count(reaction_size)?;
    let mut reactions = Vec::with_capacity(n);
    for _ in 0..n {
        let r = c.take(reaction_size)?;
        reactions.push(Reaction {
            precursor_mz: f64_at(r, 0),
            unknown_double: f64_at(r, 8),
            energy: f64_at(r, 16),
            unknown_long1: u32_at(r, 24),
            unknown_long2: u32_at(r, 28),
        });
    }
    let n = c.count(16)?;
    let mut fraction_collectors = Vec::with_capacity(n);
    for _ in 0..n {
        let r = c.take(16)?;
        fraction_collectors.push(FractionCollector {
            low_mz: f64_at(r, 0),
            high_mz: f64_at(r, 8),
        });
    }
    let n = c.count(8)?;
    let coefficients = c
        .take(n * 8)?
        .chunks_exact(8)
        .map(|b| f64_at(b, 0))
        .collect();
    let n = c.count(8)?;
    c.take(n * 8)?;
    let n = c.count(16)?;
    c.take(n * 16)?;
    if text {
        let n = c.count(2)?;
        c.take(n * 2)?;
    }
    Some(ScanEvent {
        preamble: ScanEventPreamble { bytes: preamble },
        reactions,
        fraction_collectors,
        coefficients,
    })
}

/// Reaction record sizes to try for `version`, most likely first.
fn reaction_sizes(version: u32) -> &'static [usize] {
    match version {
        0..=64 => &[32],
        65 => &[56, 48, 32],
        _ => &[56],
    }
}

/// Decode `count` events from `stream`, which must be consumed exactly.
/// Returns `None` when no known layout fits.
pub(crate) fn decode(stream: &[u8], version: u32, count: u32) -> Option<Vec<ScanEvent>> {
    let preamble_size = ScanEventPreamble::size_for_version(version);
    let text = version >= 66;
    'layout: for &reaction_size in reaction_sizes(version) {
        let mut c = Cursor {
            bytes: stream,
            pos: 0,
        };
        // A preamble is the smallest possible event.
        if (count as usize).checked_mul(preamble_size)? > stream.len() {
            return None;
        }
        let mut events = Vec::with_capacity(count as usize);
        for _ in 0..count {
            match decode_event(&mut c, preamble_size, reaction_size, text) {
                Some(event) => events.push(event),
                None => continue 'layout,
            }
        }
        if c.pos == stream.len() {
            return Some(events);
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    fn event(power: u8, reactions: &[f64], calibrators: usize, source_frag: usize) -> Vec<u8> {
        let mut b = vec![0u8; 136];
        b[6] = power;
        b[10] = u8::from(power > 1);
        b.extend((reactions.len() as u32).to_le_bytes());
        for &mz in reactions {
            let mut r = vec![0u8; 56];
            r[..8].copy_from_slice(&mz.to_le_bytes());
            r[8..16].copy_from_slice(&2f64.to_le_bytes());
            r[16..24].copy_from_slice(&30f64.to_le_bytes());
            r[24..28].copy_from_slice(&11u32.to_le_bytes());
            b.extend(r);
        }
        b.extend(1u32.to_le_bytes());
        b.extend(100f64.to_le_bytes());
        b.extend(1500f64.to_le_bytes());
        b.extend((calibrators as u32).to_le_bytes());
        for i in 0..calibrators {
            b.extend((i as f64 + 1.0).to_le_bytes());
        }
        b.extend((source_frag as u32).to_le_bytes());
        b.extend(vec![0u8; 8 * source_frag]);
        b.extend(0u32.to_le_bytes());
        let name: Vec<u16> = "PEPTIDE".encode_utf16().collect();
        b.extend((name.len() as u32).to_le_bytes());
        b.extend(name.iter().flat_map(|c| c.to_le_bytes()));
        b
    }

    #[test]
    fn decodes_variable_length_events() {
        let stream = [
            event(1, &[], 7, 0),
            event(2, &[500.0], 0, 0),
            event(2, &[500.0, 500.0], 7, 0),
            event(3, &[500.0, 300.0], 7, 1),
            event(2, &[400.0, 410.0, 420.0, 430.0], 5, 1),
        ]
        .concat();
        let events = decode(&stream, 66, 5).unwrap();
        assert_eq!(events.len(), 5);
        assert_eq!(events[0].coefficients.len(), 7);
        assert!(events[1].coefficients.is_empty());
        assert_eq!(events[2].reactions.len(), 2);
        assert_eq!(events[3].preamble.bytes[6], 3);
        assert_eq!(events[3].reactions[1].precursor_mz, 300.0);
        assert_eq!(events[4].reactions[3].precursor_mz, 430.0);
        assert_eq!(events[4].fraction_collectors[0].high_mz, 1500.0);
    }

    #[test]
    fn rejects_streams_that_do_not_end_exactly() {
        let stream = [event(1, &[], 7, 0), event(2, &[500.0], 7, 0)].concat();
        assert!(decode(&stream[..stream.len() - 1], 66, 2).is_none());
        let mut longer = stream.clone();
        longer.push(0);
        assert!(decode(&longer, 66, 2).is_none());
        assert!(decode(&stream, 66, 3).is_none());
        let mut bad_count = stream;
        bad_count[136..140].copy_from_slice(&u32::MAX.to_le_bytes());
        assert!(decode(&bad_count, 66, 2).is_none());
    }
}
