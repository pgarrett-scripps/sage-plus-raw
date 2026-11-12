use crate::binary_reader::BinaryReader;
use crate::error::{Error, Result};
use std::io::{Read, Seek, SeekFrom};

/// A single centroid peak (m/z + abundance).
#[derive(Debug, Clone)]
pub struct Peak {
    pub mz: f64,
    pub abundance: f32,
}

/// A node in the per-scan noise-vs-m/z function carried by FT scans.
///
/// Thermo stores noise and baseline as a piecewise-linear function of m/z
/// (a few dozen nodes), not per peak. Per-peak noise/baseline are recovered
/// by interpolating this function at each peak's m/z (see
/// [`ScanDataPacket::noise_at`]).
#[derive(Debug, Clone, Copy)]
pub struct NoiseNode {
    pub mz: f32,
    pub noise: f32,
    pub baseline: f32,
}

/// A contiguous chunk of profile signal data.
#[derive(Debug)]
pub struct ProfileChunk {
    pub first_bin: u32,
    pub signal: Vec<f32>,
    pub fudge: Option<f32>,
}

/// Profile spectrum data.
#[derive(Debug)]
pub struct Profile {
    pub first_value: f64,
    pub step: f64,
    pub peak_count: u32,
    pub nbins: u32,
    pub chunks: Vec<ProfileChunk>,
}

/// The header of a scan data packet (40 bytes).
#[derive(Debug)]
pub struct PacketHeader {
    pub profile_size: u32,
    pub peak_list_size: u32,
    pub layout: u32,
    pub descriptor_list_size: u32,
    pub unknown_stream_size: u32,
    pub triplet_stream_size: u32,
    pub low_mz: f32,
    pub high_mz: f32,
}

/// A complete scan data packet.
#[derive(Debug)]
pub struct ScanDataPacket {
    pub header: PacketHeader,
    pub profile: Option<Profile>,
    pub peaks: Vec<Peak>,
    /// Per-peak resolution, aligned 1:1 with `peaks`. Empty when the scan
    /// carries no FT label data (e.g. ion-trap scans) or its layout did not
    /// match the expected encoding.
    pub resolutions: Vec<f32>,
    /// Nodes of the scan's noise-vs-m/z function. Empty when absent. Use
    /// [`Self::noise_at`] to evaluate per-peak noise/baseline.
    pub noise_nodes: Vec<NoiseNode>,
}

impl ScanDataPacket {
    pub(crate) fn read<R: Read + Seek>(
        r: &mut BinaryReader<R>,
        scratch: &mut Vec<u8>,
    ) -> Result<Self> {
        let header = PacketHeader::read(r)?;

        // Profile data
        let profile = if header.profile_size > 0 {
            Some(Profile::read(r, header.layout, scratch)?)
        } else {
            None
        };

        let peaks = Self::read_peaks(r, &header, scratch)?;
        let (resolutions, noise_nodes) = Self::read_labels(r, &header, peaks.len(), scratch)?;

        Ok(Self {
            header,
            profile,
            peaks,
            resolutions,
            noise_nodes,
        })
    }

    /// Read peaks and FT label data while skipping the (potentially large)
    /// profile signal. This is the fast path when only centroids and their
    /// labels (resolution / noise / baseline) are needed.
    pub(crate) fn read_skip_profile<R: Read + Seek>(
        r: &mut BinaryReader<R>,
        scratch: &mut Vec<u8>,
    ) -> Result<Self> {
        let header = PacketHeader::read(r)?;

        if header.profile_size > 0 {
            r.skip((header.profile_size as usize) * 4)?;
        }

        let peaks = Self::read_peaks(r, &header, scratch)?;
        let (resolutions, noise_nodes) = Self::read_labels(r, &header, peaks.len(), scratch)?;

        Ok(Self {
            header,
            profile: None,
            peaks,
            resolutions,
            noise_nodes,
        })
    }

    /// Read only the centroided peak list, skipping the profile data and the
    /// trailing label streams. This is 2-10× faster than [`Self::read`] for
    /// high-resolution Orbitrap scans where profile_size can be tens of
    /// thousands of 4-byte words.
    pub(crate) fn read_peaks_only<R: Read + Seek>(
        r: &mut BinaryReader<R>,
        scratch: &mut Vec<u8>,
    ) -> Result<Vec<Peak>> {
        let (block, wide_mz) = Self::read_peaks_only_block(r, scratch)?;
        Ok(peaks_from_block(block, wide_mz))
    }

    /// [`Self::read_peaks_only`] as the undecoded peak-list block, for a
    /// caller that fills its own arrays from it ([`for_each_peak`]) instead
    /// of taking a `Vec<Peak>` it would only copy out of.
    pub(crate) fn read_peaks_only_block<'a, R: Read + Seek>(
        r: &mut BinaryReader<R>,
        scratch: &'a mut Vec<u8>,
    ) -> Result<(&'a [u8], bool)> {
        let header = PacketHeader::read(r)?;

        // Skip profile data instead of decoding it.
        if header.profile_size > 0 {
            r.skip((header.profile_size as usize) * 4)?;
        }

        Self::read_peak_block(r, &header, scratch)
    }

    /// Decode the centroid peak list. Assumes the cursor sits at the start of
    /// the peak list (i.e. profile data has already been read or skipped).
    fn read_peaks<R: Read + Seek>(
        r: &mut BinaryReader<R>,
        header: &PacketHeader,
        scratch: &mut Vec<u8>,
    ) -> Result<Vec<Peak>> {
        let (block, wide_mz) = Self::read_peak_block(r, header, scratch)?;
        Ok(peaks_from_block(block, wide_mz))
    }

    /// Read the scan's peak list as bytes, saying whether its m/z fields are
    /// f64 (`true`) or f32. Assumes the cursor sits at the start of the peak
    /// list (i.e. profile data has already been read or skipped).
    fn read_peak_block<'a, R: Read + Seek>(
        r: &mut BinaryReader<R>,
        header: &PacketHeader,
        scratch: &'a mut Vec<u8>,
    ) -> Result<(&'a [u8], bool)> {
        // Layout bit 16 (0x10000) means m/z is f64 instead of f32.
        let wide_mz = header.layout & 0x10000 != 0;
        if header.peak_list_size == 0 {
            return Ok((&[], wide_mz));
        }
        let count = r.read_u32()? as usize;
        // The whole peak list in one read, then parsed from the bytes: a
        // field at a time is two reads per peak, around 2,400 per scan.
        let block = r.read_items_into(scratch, count, peak_item_size(wide_mz))?;
        Ok((block, wide_mz))
    }

    /// Decode the per-peak label streams that follow the centroid peak list:
    /// descriptor (peak index, skipped), unknown (per-peak resolution), and
    /// triplet (noise-vs-m/z function). Assumes the cursor sits immediately
    /// after the peak list.
    fn read_labels<R: Read + Seek>(
        r: &mut BinaryReader<R>,
        header: &PacketHeader,
        n_peaks: usize,
        scratch: &mut Vec<u8>,
    ) -> Result<(Vec<f32>, Vec<NoiseNode>)> {
        // Descriptor stream: one u32 per peak (a peak index), not label data.
        //
        // Cast to usize *before* multiplying: `descriptor_list_size` is an
        // untrusted u32 read straight from the file, and multiplying by 4 in
        // u32 arithmetic first can overflow (and panic, in a debug/fuzz
        // build) for values above ~1.07 billion. Widening first is safe
        // because `usize` is 64-bit on every platform this crate targets.
        if header.descriptor_list_size > 0 {
            r.skip(header.descriptor_list_size as usize * 4)?;
        }

        // Unknown stream: for FT scans this is `[count_word, resolution...]`
        // with one f32 resolution per centroid peak. Decode it only when the
        // size matches that layout; otherwise skip it so other encodings are
        // left untouched.
        let resolutions = if n_peaks > 0 && header.unknown_stream_size as usize == n_peaks + 1 {
            let _count = r.read_u32()?;
            let block = r.read_items_into(scratch, n_peaks, 4)?;
            let mut res = Vec::with_capacity(n_peaks);
            for item in block.chunks_exact(4) {
                res.push(f32::from_le_bytes([item[0], item[1], item[2], item[3]]));
            }
            res
        } else {
            // Same overflow hazard as the descriptor stream above: widen
            // before multiplying.
            if header.unknown_stream_size > 0 {
                r.skip(header.unknown_stream_size as usize * 4)?;
            }
            Vec::new()
        };

        // Triplet stream: the noise-vs-m/z function as (m/z, noise, baseline)
        // f32 nodes. Per-peak noise/baseline come from interpolating this at
        // the peak m/z (see `noise_at`).
        let node_count = header.triplet_stream_size / 3;
        let block = r.read_items_into(scratch, node_count as usize, 12)?;
        let mut noise_nodes = Vec::with_capacity(node_count as usize);
        for item in block.chunks_exact(12) {
            noise_nodes.push(NoiseNode {
                mz: f32::from_le_bytes([item[0], item[1], item[2], item[3]]),
                noise: f32::from_le_bytes([item[4], item[5], item[6], item[7]]),
                baseline: f32::from_le_bytes([item[8], item[9], item[10], item[11]]),
            });
        }
        // Skip any trailing words if the stream size is not a multiple of 3.
        // The remainder here is always < 3 (by construction of `node_count`
        // as a floor division), so widening before multiplying isn't
        // strictly required for overflow safety, but it matches the pattern
        // used everywhere else in this function for consistency.
        let consumed = node_count * 3;
        if header.triplet_stream_size > consumed {
            r.skip((header.triplet_stream_size - consumed) as usize * 4)?;
        }

        Ok((resolutions, noise_nodes))
    }

    /// Linearly interpolate `(noise, baseline)` at `mz` from the scan's
    /// noise-vs-m/z function. Returns `None` when the scan carries no noise
    /// nodes. Queries outside the node range clamp to the nearest endpoint.
    pub fn noise_at(&self, mz: f64) -> Option<(f32, f32)> {
        let nodes = &self.noise_nodes;
        if nodes.is_empty() {
            return None;
        }
        let x = mz as f32;
        if x <= nodes[0].mz {
            return Some((nodes[0].noise, nodes[0].baseline));
        }
        for w in nodes.windows(2) {
            let (lo, hi) = (&w[0], &w[1]);
            if x <= hi.mz {
                let span = hi.mz - lo.mz;
                let f = if span > 0.0 { (x - lo.mz) / span } else { 0.0 };
                let noise = lo.noise + f * (hi.noise - lo.noise);
                let baseline = lo.baseline + f * (hi.baseline - lo.baseline);
                return Some((noise, baseline));
            }
        }
        let last = nodes[nodes.len() - 1];
        Some((last.noise, last.baseline))
    }
}

/// Size on disk of one peak-list item: an f32 intensity and an m/z that is
/// f64 in the wide layout and f32 otherwise.
pub(crate) fn peak_item_size(wide_mz: bool) -> usize {
    if wide_mz {
        12
    } else {
        8
    }
}

/// Parse every peak of a peak-list block, handing `(mz, intensity)` to
/// `visit` in file order.
///
/// This is the one place peak-list bytes become numbers, so a caller filling
/// its own m/z and intensity arrays gets exactly the values
/// [`ScanDataPacket::read_peaks_only`] puts in its `Vec<Peak>`.
pub(crate) fn for_each_peak(block: &[u8], wide_mz: bool, mut visit: impl FnMut(f64, f32)) {
    // Each item is read at fixed offsets in a chunk of known length, so
    // the bounds checks fold away and the loop is one load per field.
    if wide_mz {
        for item in block.chunks_exact(12) {
            visit(
                f64::from_le_bytes([
                    item[0], item[1], item[2], item[3], item[4], item[5], item[6], item[7],
                ]),
                f32::from_le_bytes([item[8], item[9], item[10], item[11]]),
            );
        }
    } else {
        for item in block.chunks_exact(8) {
            visit(
                f64::from(f32::from_le_bytes([item[0], item[1], item[2], item[3]])),
                f32::from_le_bytes([item[4], item[5], item[6], item[7]]),
            );
        }
    }
}

/// Collect a peak-list block into a `Vec<Peak>`.
fn peaks_from_block(block: &[u8], wide_mz: bool) -> Vec<Peak> {
    let mut peaks = Vec::with_capacity(block.len() / peak_item_size(wide_mz));
    for_each_peak(block, wide_mz, |mz, abundance| {
        peaks.push(Peak { mz, abundance })
    });
    peaks
}

impl PacketHeader {
    fn read<R: Read + Seek>(r: &mut BinaryReader<R>) -> Result<Self> {
        let _unk1 = r.read_u32()?;
        let profile_size = r.read_u32()?;
        let peak_list_size = r.read_u32()?;
        let layout = r.read_u32()?;
        let descriptor_list_size = r.read_u32()?;
        let unknown_stream_size = r.read_u32()?;
        let triplet_stream_size = r.read_u32()?;
        let _unk2 = r.read_u32()?;
        let low_mz = r.read_f32()?;
        let high_mz = r.read_f32()?;

        Ok(Self {
            profile_size,
            peak_list_size,
            layout,
            descriptor_list_size,
            unknown_stream_size,
            triplet_stream_size,
            low_mz,
            high_mz,
        })
    }
}

impl Profile {
    fn read<R: Read + Seek>(
        r: &mut BinaryReader<R>,
        layout: u32,
        scratch: &mut Vec<u8>,
    ) -> Result<Self> {
        let first_value = r.read_f64()?;
        let step = r.read_f64()?;
        let peak_count = r.read_u32()?;
        let nbins = r.read_u32()?;

        let has_fudge = layout & 0xFF != 0;
        // Minimum on-disk size of one chunk: first_bin(4) + chunk_nbins(4).
        r.check_count(peak_count as u64, 8)?;
        let mut chunks = Vec::with_capacity(peak_count as usize);
        for _ in 0..peak_count {
            let first_bin = r.read_u32()?;
            let chunk_nbins = r.read_u32()?;
            let fudge = if has_fudge { Some(r.read_f32()?) } else { None };
            let block = r.read_items_into(scratch, chunk_nbins as usize, 4)?;
            let mut signal = Vec::with_capacity(chunk_nbins as usize);
            for item in block.chunks_exact(4) {
                signal.push(f32::from_le_bytes([item[0], item[1], item[2], item[3]]));
            }
            chunks.push(ProfileChunk {
                first_bin,
                signal,
                fudge,
            });
        }

        Ok(Self {
            first_value,
            step,
            peak_count,
            nbins,
            chunks,
        })
    }
}

impl Profile {
    /// Convert profile bins to (mz, intensity) pairs using the conversion coefficients.
    ///
    /// When the coefficient count is not a known calibration layout (see
    /// [`freq_to_mz`]), every m/z is `NaN`.
    pub fn to_mz_intensity(&self, coefficients: &[f64]) -> Vec<(f64, f64)> {
        // Size the allocation from the actually-parsed (and therefore already
        // bounded) chunk signal lengths, not from `self.nbins`: that field is
        // an untrusted file-provided value not otherwise validated against
        // the real signal data, so using it directly here would reopen the
        // same unbounded-allocation issue the read-path checks close off.
        let mut result = Vec::with_capacity(self.signal_len());
        result.extend(self.mz_intensity_iter(coefficients));
        result
    }

    /// Total parsed signal length across chunks: how many pairs
    /// [`Self::mz_intensity_iter`] yields.
    pub(crate) fn signal_len(&self) -> usize {
        self.chunks.iter().map(|c| c.signal.len()).sum()
    }

    /// The same (m/z, intensity) pairs as [`Self::to_mz_intensity`], without
    /// collecting them, for callers filling their own buffers.
    pub(crate) fn mz_intensity_iter<'a>(
        &'a self,
        coefficients: &'a [f64],
    ) -> impl Iterator<Item = (f64, f64)> + 'a {
        self.chunks.iter().flat_map(move |chunk| {
            chunk.signal.iter().enumerate().map(move |(i, &intensity)| {
                let bin_global = chunk.first_bin as f64 + i as f64;
                let freq = self.first_value + bin_global * self.step;
                // The per-chunk fudge is an additive m/z correction applied
                // after frequency-to-m/z conversion, not a frequency offset.
                // Applying it in the frequency domain shifts Orbitrap profile
                // bins by ~1e-3 m/z relative to the centroid list (verified
                // against centroids across Q Exactive / Exploris / Eclipse).
                let mz = freq_to_mz(freq, coefficients) + chunk.fudge.unwrap_or(0.0) as f64;
                (mz, intensity as f64)
            })
        })
    }
}

/// Convert frequency to m/z using the conversion coefficients from the scan event.
///
/// The coefficient array includes metadata prefix values:
/// - nparam=4 (LTQ-FT/ICR): [unknown, A, B, C] → Mz = A + B/f + C/f²
/// - nparam=5 (Orbitrap v66): [unk0, unk1, A, B, C] → Mz = A + B/f² + C/f⁴
/// - nparam=7 (Orbitrap): [unknown, I, A, B, C, D, E] → Mz = A + B/f² + C/f⁴
///
/// Any other non-zero coefficient count is an unknown calibration layout and
/// returns `f64::NAN`, so an unconverted frequency is never mistaken for m/z.
/// Callers must treat non-finite results as "no m/z available".
pub fn freq_to_mz(freq: f64, coefficients: &[f64]) -> f64 {
    if freq == 0.0 {
        return 0.0;
    }
    match coefficients.len() {
        0 => freq, // No conversion (already m/z domain, e.g. ITMS)
        4 => {
            // LTQ-FT / ICR: Mz = A + B/f + C/f²
            let (a, b, c) = (coefficients[1], coefficients[2], coefficients[3]);
            a + b / freq + c / (freq * freq)
        }
        5 => {
            // Orbitrap v66: Mz = A + B/f² + C/f⁴
            let (a, b, c) = (coefficients[2], coefficients[3], coefficients[4]);
            let f2 = freq * freq;
            a + b / f2 + c / (f2 * f2)
        }
        7 => {
            // Orbitrap: Mz = A + B/f² + C/f⁴
            let (a, b, c) = (coefficients[2], coefficients[3], coefficients[4]);
            let f2 = freq * freq;
            a + b / f2 + c / (f2 * f2)
        }
        _ => f64::NAN,
    }
}

/// Read a flat-peak scan (TSQ/SRM format).
///
/// In this format, the scan data stream at `data_addr` contains variable-length records.
/// Each scan index entry's `offset` field holds the **cumulative end byte offset** within
/// the data stream. Peaks are stored as contiguous (f32 mz, f32 intensity) pairs at the
/// end of each record, followed by `peak_count` flag bytes (1 byte per peak).
///
/// `peak_count` is typically `data_size - 1`.
pub fn read_flat_peaks<R: Read + Seek>(
    source: &mut R,
    data_addr: u64,
    cum_end: u64,
    data_size: u32,
) -> Result<Vec<Peak>> {
    if data_size <= 1 {
        return Ok(Vec::new());
    }

    // Try peak_count = data_size - 1 first, then data_size - 2 as fallback.
    // Each peak occupies 9 bytes total: 8 bytes (f32 mz + f32 int) + 1 flag byte.
    // The peaks section is at the end of the record.
    let mut scratch = Vec::new();
    for subtract in [1u32, 2] {
        if data_size <= subtract {
            continue;
        }
        let peak_count = (data_size - subtract) as usize;
        let peak_section_bytes = peak_count as u64 * 9;
        if peak_section_bytes > cum_end {
            continue;
        }
        // `data_addr` and `cum_end` are both untrusted u64 values read
        // straight from the file (a RawFileInfo address and a ScanIndexEntry
        // offset respectively); saturate instead of panicking on overflow
        // for a corrupt/crafted combination of the two.
        let peaks_start = data_addr
            .saturating_add(cum_end)
            .saturating_sub(peak_section_bytes);
        source.seek(SeekFrom::Start(peaks_start))?;
        let mut r = BinaryReader::new(&mut *source);

        // The whole peak section in one read, then parsed from the bytes,
        // as the PacketHeader peak list is: a field at a time is two reads
        // per peak. `read_items_into` applies the same `check_count`
        // allocation bound (`peak_count` items of 8 bytes) this loop did.
        let block = r.read_items_into(&mut scratch, peak_count, 8)?;
        let mut peaks = Vec::with_capacity(peak_count);
        for item in block.chunks_exact(8) {
            peaks.push(Peak {
                mz: f64::from(f32::from_le_bytes([item[0], item[1], item[2], item[3]])),
                abundance: f32::from_le_bytes([item[4], item[5], item[6], item[7]]),
            });
        }

        // Validate: first peak mz should be a plausible mass value (or zero for empty transitions)
        let looks_valid = if let Some(first) = peaks.first() {
            first.mz == 0.0 || (first.mz > 10.0 && first.mz < 10_000.0)
        } else {
            true
        };

        if looks_valid {
            return Ok(peaks);
        }
    }

    // Neither worked; return empty
    Err(Error::UnexpectedEof {
        offset: data_addr.saturating_add(cum_end),
        needed: 0,
    })
}

/// Read a flat-peak scan in v66 (TSQ Quantiva / TSQ Altis) SRM format.
///
/// In this format, the scan data stream contains fixed-size records.
/// Each scan index entry's `offset` field holds the **start byte offset** within
/// the stream (not cumulative end), and `record_size` is the number of bytes per record.
///
/// Record layout:
/// - bytes 0-3: u32 `n_peaks` (number of active SRM transitions in this window)
/// - bytes 4-31: other header fields (skipped)
/// - bytes 32..32+n_peaks*8: m/z window table, one (lo_mz: f32, hi_mz: f32) pair per peak
/// - bytes 32+n_peaks*8..: peak data, one (channel_idx: u32, mz: f32, intensity: f32) per peak
pub fn read_scan_srm_v66<R: Read + Seek>(
    source: &mut R,
    data_addr: u64,
    start_offset: u64,
    _record_size: u32,
) -> Result<Vec<Peak>> {
    // Both are untrusted u64 values read from the file; saturate
    // instead of panicking on overflow for a corrupt combination.
    let abs_start = data_addr.saturating_add(start_offset);
    source.seek(SeekFrom::Start(abs_start))?;
    let mut r = BinaryReader::new(source);

    // n_peaks at byte 0
    let n_peaks = r.read_u32()? as usize;
    if n_peaks == 0 {
        return Ok(Vec::new());
    }

    // Skip remaining header: bytes 4-31 (28 bytes)
    r.skip(28)?;

    // Skip m/z window table: n_peaks × 8 bytes (lo_mz f32 + hi_mz f32 per channel)
    r.skip(n_peaks * 8)?;

    // Read peak records: (u32 channel_idx, f32 mz, f32 intensity) × n_peaks.
    // The whole table in one read, then parsed from the bytes;
    // `read_items_into` applies the same `check_count` allocation bound
    // (`n_peaks` items of 12 bytes) the field-at-a-time loop did.
    let mut scratch = Vec::new();
    let block = r.read_items_into(&mut scratch, n_peaks, 12)?;
    let mut peaks = Vec::with_capacity(n_peaks);
    for item in block.chunks_exact(12) {
        // Bytes 0-3 are the channel index, which this reader does not use.
        peaks.push(Peak {
            mz: f64::from(f32::from_le_bytes([item[4], item[5], item[6], item[7]])),
            abundance: f32::from_le_bytes([item[8], item[9], item[10], item[11]]),
        });
    }

    Ok(peaks)
}

/// Read the Q3 isolation window table from an SRM v66 scan record.
///
/// Returns one `(lo_mz, hi_mz)` pair per active transition channel,
/// in channel order.
pub fn read_scan_srm_v66_windows<R: Read + Seek>(
    source: &mut R,
    data_addr: u64,
    start_offset: u64,
) -> Result<Vec<(f32, f32)>> {
    // Both are untrusted u64 values read from the file; saturate
    // instead of panicking on overflow for a corrupt combination.
    let abs_start = data_addr.saturating_add(start_offset);
    source.seek(SeekFrom::Start(abs_start))?;
    let mut r = BinaryReader::new(source);

    // n_peaks at byte 0
    let n_peaks = r.read_u32()? as usize;
    if n_peaks == 0 {
        return Ok(Vec::new());
    }

    // Skip remaining header: bytes 4-31 (28 bytes)
    r.skip(28)?;

    // Read m/z window table: n_peaks × 8 bytes (lo_mz f32, hi_mz f32).
    // One read for the table, then parsed from the bytes;
    // `read_items_into` applies the same `check_count` allocation bound
    // (`n_peaks` items of 8 bytes) the field-at-a-time loop did.
    let mut scratch = Vec::new();
    let block = r.read_items_into(&mut scratch, n_peaks, 8)?;
    let mut windows = Vec::with_capacity(n_peaks);
    for item in block.chunks_exact(8) {
        windows.push((
            f32::from_le_bytes([item[0], item[1], item[2], item[3]]),
            f32::from_le_bytes([item[4], item[5], item[6], item[7]]),
        ));
    }

    Ok(windows)
}

/// Search the pre-data method/transition table for a v63 SRM transition record
/// matching the given Q3 center mass.
///
/// v63 (TSQ Quantum/Vantage) transition table layout: 72 bytes per channel record.
/// Relevant fields (all f64 little-endian):
///   - [+ 0] active-channel flag (1.0 for first channel of each precursor)
///   - [+ 8] unknown
///   - [+16] Q1 precursor mass (m/z)   ← returned
///   - [+24] Q3 center mass (m/z)      ← anchor for search
///   - [+32] Q3 window width (Da)      ← returned
///   - [+40] dwell time (s)
///   - [+48] collision energy (eV)     ← returned
///
/// Returns `(Q1, Q3_width, CE_eV)` if a plausible match is found.
pub fn search_v63_transition(data: &[u8], q3_center_target: f64) -> Option<(f64, f64, f64)> {
    let end = data.len().saturating_sub(32);
    for j in 8..end {
        if j + 8 > data.len() {
            break;
        }
        let v = f64::from_le_bytes(data[j..j + 8].try_into().ok()?);
        if (v - q3_center_target).abs() > 0.002 {
            continue;
        }
        // Candidate Q3_center at position j. Q1 is 8 bytes before.
        let q1 = f64::from_le_bytes(data[j - 8..j].try_into().ok()?);
        if !q1.is_finite() || !(50.0..=3000.0).contains(&q1) {
            continue;
        }
        // Q3_width is 8 bytes after Q3_center.
        if j + 16 > data.len() {
            continue;
        }
        let q3w = f64::from_le_bytes(data[j + 8..j + 16].try_into().ok()?);
        if !q3w.is_finite() || !(0.01..=10.0).contains(&q3w) {
            continue;
        }
        // CE is 24 bytes after Q3_center.
        if j + 32 > data.len() {
            continue;
        }
        let ce = f64::from_le_bytes(data[j + 24..j + 32].try_into().ok()?);
        if !ce.is_finite() || !(0.1..=300.0).contains(&ce) {
            continue;
        }
        return Some((q1, q3w, ce));
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn freq_to_mz_passes_through_without_coefficients() {
        assert_eq!(freq_to_mz(500.0, &[]), 500.0);
        assert_eq!(freq_to_mz(0.0, &[1.0; 3]), 0.0);
    }

    #[test]
    fn freq_to_mz_is_nan_for_unknown_coefficient_counts() {
        for n in [1, 2, 3, 6, 8, 12] {
            let coeffs = vec![1.0; n];
            assert!(freq_to_mz(123_456.0, &coeffs).is_nan(), "nparam={n}");
        }
        // Known layouts stay finite.
        for n in [4, 5, 7] {
            let coeffs = vec![1.0; n];
            assert!(freq_to_mz(123_456.0, &coeffs).is_finite(), "nparam={n}");
        }
    }

    fn header_with(
        descriptor_list_size: u32,
        unknown_stream_size: u32,
        triplet_stream_size: u32,
    ) -> PacketHeader {
        PacketHeader {
            profile_size: 0,
            peak_list_size: 0,
            layout: 0,
            descriptor_list_size,
            unknown_stream_size,
            triplet_stream_size,
            low_mz: 0.0,
            high_mz: 0.0,
        }
    }

    /// Regression test for a fuzzer-found crash: `descriptor_list_size` and
    /// `unknown_stream_size` are untrusted u32 values read straight from the
    /// file. Multiplying by 4 in u32 arithmetic (rather than widening first)
    /// overflows - and panics under debug assertions - for values at or
    /// above `u32::MAX / 4 + 1`. `read_labels` must reject or skip past
    /// these cleanly instead of panicking.
    #[test]
    fn huge_descriptor_list_size_does_not_overflow() {
        let header = header_with(u32::MAX / 4 + 1, 0, 0);
        let mut r = BinaryReader::new(Cursor::new(vec![0u8; 8]));
        // Skipping past EOF is allowed (only a subsequent read would fail);
        // the point of this test is that it must not panic.
        let _ = ScanDataPacket::read_labels(&mut r, &header, 0, &mut Vec::new());
    }

    #[test]
    fn huge_unknown_stream_size_does_not_overflow() {
        // n_peaks = 0 so the `unknown_stream_size == n_peaks + 1` decode
        // branch is not taken and the plain skip path (the one that
        // crashed) runs instead.
        let header = header_with(0, u32::MAX / 4 + 1, 0);
        let mut r = BinaryReader::new(Cursor::new(vec![0u8; 8]));
        let _ = ScanDataPacket::read_labels(&mut r, &header, 0, &mut Vec::new());
    }

    #[test]
    fn read_labels_with_no_streams_returns_empty() {
        let header = header_with(0, 0, 0);
        let mut r = BinaryReader::new(Cursor::new(Vec::new()));
        let (resolutions, noise_nodes) =
            ScanDataPacket::read_labels(&mut r, &header, 0, &mut Vec::new()).unwrap();
        assert!(resolutions.is_empty());
        assert!(noise_nodes.is_empty());
    }

    #[test]
    fn read_labels_decodes_resolutions_when_size_matches() {
        // unknown_stream_size == n_peaks + 1 selects the resolution-decoding
        // branch: a leading count word followed by one f32 per peak.
        let header = header_with(0, 3, 0);
        let mut bytes = 2u32.to_le_bytes().to_vec(); // count word
        bytes.extend_from_slice(&1.5f32.to_le_bytes());
        bytes.extend_from_slice(&2.5f32.to_le_bytes());
        let mut r = BinaryReader::new(Cursor::new(bytes));
        let (resolutions, _) =
            ScanDataPacket::read_labels(&mut r, &header, 2, &mut Vec::new()).unwrap();
        assert_eq!(resolutions, vec![1.5, 2.5]);
    }

    #[test]
    fn read_labels_decodes_noise_nodes() {
        let header = header_with(0, 0, 3); // one (mz, noise, baseline) node
        let mut bytes = 100.0f32.to_le_bytes().to_vec();
        bytes.extend_from_slice(&5.0f32.to_le_bytes());
        bytes.extend_from_slice(&1.0f32.to_le_bytes());
        let mut r = BinaryReader::new(Cursor::new(bytes));
        let (_, noise_nodes) =
            ScanDataPacket::read_labels(&mut r, &header, 0, &mut Vec::new()).unwrap();
        assert_eq!(noise_nodes.len(), 1);
        assert_eq!(noise_nodes[0].mz, 100.0);
        assert_eq!(noise_nodes[0].noise, 5.0);
        assert_eq!(noise_nodes[0].baseline, 1.0);
    }

    #[test]
    fn read_peaks_rejects_implausible_count_before_allocating() {
        // peak_list_size must be nonzero or read_peaks short-circuits to an
        // empty Vec before ever reading a count.
        let mut header = header_with(0, 0, 0);
        header.peak_list_size = 1;
        let mut bytes = u32::MAX.to_le_bytes().to_vec(); // declared peak count
        bytes.extend_from_slice(&[0u8; 4]); // far short of what's needed
        let mut r = BinaryReader::new(Cursor::new(bytes));
        let err = ScanDataPacket::read_peaks(&mut r, &header, &mut Vec::new()).unwrap_err();
        assert!(matches!(err, Error::AllocationTooLarge { .. }));
    }

    /// Peak-list bytes: a count word followed by `count` items, narrow
    /// (f32 m/z) or wide (f64 m/z) as the layout bit says.
    fn peak_list_bytes(peaks: &[(f64, f32)], wide_mz: bool) -> Vec<u8> {
        let mut bytes = (peaks.len() as u32).to_le_bytes().to_vec();
        for (mz, abundance) in peaks {
            if wide_mz {
                bytes.extend_from_slice(&mz.to_le_bytes());
            } else {
                bytes.extend_from_slice(&(*mz as f32).to_le_bytes());
            }
            bytes.extend_from_slice(&abundance.to_le_bytes());
        }
        bytes
    }

    #[test]
    fn read_peaks_decodes_the_block_the_same_in_both_layouts() {
        for wide_mz in [false, true] {
            let mut header = header_with(0, 0, 0);
            header.peak_list_size = 1;
            header.layout = if wide_mz { 0x10000 } else { 0 };
            let expected = [(500.5f64, 10.0f32), (700.25, 2.5), (1000.125, 0.0)];
            let bytes = peak_list_bytes(&expected, wide_mz);
            let item_size = if wide_mz { 12 } else { 8 };
            assert_eq!(bytes.len(), 4 + expected.len() * item_size);
            let mut r = BinaryReader::new(Cursor::new(bytes));
            let peaks = ScanDataPacket::read_peaks(&mut r, &header, &mut Vec::new()).unwrap();
            assert_eq!(peaks.len(), expected.len());
            for (peak, (mz, abundance)) in peaks.iter().zip(expected) {
                // Narrow m/z is an f32 widened, which is what the old
                // field-at-a-time read produced too.
                let mz = if wide_mz { mz } else { f64::from(mz as f32) };
                assert_eq!((peak.mz, peak.abundance), (mz, abundance));
            }
        }
    }

    #[test]
    fn read_peaks_truncated_reports_the_peak_it_could_not_finish() {
        // Two whole peaks and half of a third, with the count claiming
        // three: the error is the one a peak-at-a-time read gave, at the
        // offset of the peak that could not be read whole.
        let mut header = header_with(0, 0, 0);
        header.peak_list_size = 1;
        let mut bytes = peak_list_bytes(&[(500.0, 1.0), (600.0, 2.0)], false);
        bytes[0] = 3;
        bytes.extend_from_slice(&700.0f32.to_le_bytes());
        // Told the file is longer than this packet, as a per-scan reader is,
        // so the allocation cap passes and the short read is what fails.
        let mut r = BinaryReader::at(Cursor::new(bytes), 0, 64);
        let error = ScanDataPacket::read_peaks(&mut r, &header, &mut Vec::new()).unwrap_err();
        assert!(
            matches!(
                error,
                Error::UnexpectedEof {
                    offset: 20,
                    needed: 8
                }
            ),
            "got {error:?}"
        );
    }

    #[test]
    fn read_peaks_empty_when_peak_list_size_zero() {
        let header = header_with(0, 0, 0);
        let mut r = BinaryReader::new(Cursor::new(Vec::new()));
        assert_eq!(
            ScanDataPacket::read_peaks(&mut r, &header, &mut Vec::new())
                .unwrap()
                .len(),
            0
        );
    }

    #[test]
    fn read_flat_peaks_decodes_valid_record() {
        // data_addr=0, cum_end = end offset of this scan's record. One peak
        // (mz=500.0, abundance=10.0) plus one trailing flag byte, so
        // data_size = peak_count(1) + 1 = 2.
        let mut bytes = vec![0u8; 100]; // leading padding
        let peak_start = bytes.len();
        bytes.extend_from_slice(&500.0f32.to_le_bytes());
        bytes.extend_from_slice(&10.0f32.to_le_bytes());
        bytes.push(0); // flag byte accounted for by `data_size - 1`
        let cum_end = (bytes.len()) as u64;
        let mut cursor = Cursor::new(bytes);
        let peaks = read_flat_peaks(&mut cursor, 0, cum_end, 2).unwrap();
        assert_eq!(peaks.len(), 1);
        assert_eq!(peaks[0].mz, 500.0);
        assert_eq!(peaks[0].abundance, 10.0);
        let _ = peak_start; // silence unused-var warning if layout changes
    }

    #[test]
    fn read_flat_peaks_rejects_implausible_peak_count() {
        // data_size implies a peak_count that could never fit in the file.
        let bytes = vec![0u8; 16];
        let mut cursor = Cursor::new(bytes);
        let err = read_flat_peaks(&mut cursor, 0, 16, u32::MAX).unwrap_err();
        // Either our new allocation guard or the pre-existing
        // "peak_section_bytes > cum_end" guard is an acceptable rejection
        // reason; a panic or successful huge allocation is not.
        match err {
            Error::AllocationTooLarge { .. } | Error::UnexpectedEof { .. } => {}
            other => panic!("unexpected error variant: {other:?}"),
        }
    }

    #[test]
    fn read_scan_srm_v66_decodes_valid_record() {
        // n_peaks=1, 28 bytes of skipped header, one (lo,hi) window pair,
        // then one (channel, mz, intensity) peak record.
        let mut bytes = 1u32.to_le_bytes().to_vec();
        bytes.extend_from_slice(&[0u8; 28]);
        bytes.extend_from_slice(&100.0f32.to_le_bytes()); // window lo
        bytes.extend_from_slice(&110.0f32.to_le_bytes()); // window hi
        bytes.extend_from_slice(&7u32.to_le_bytes()); // channel idx
        bytes.extend_from_slice(&105.0f32.to_le_bytes()); // mz
        bytes.extend_from_slice(&42.0f32.to_le_bytes()); // intensity
        let mut cursor = Cursor::new(bytes);
        let peaks = read_scan_srm_v66(&mut cursor, 0, 0, 0).unwrap();
        assert_eq!(peaks.len(), 1);
        assert_eq!(peaks[0].mz, 105.0);
        assert_eq!(peaks[0].abundance, 42.0);
    }

    #[test]
    fn read_scan_srm_v66_zero_peaks_is_empty() {
        let bytes = 0u32.to_le_bytes().to_vec();
        let mut cursor = Cursor::new(bytes);
        assert!(read_scan_srm_v66(&mut cursor, 0, 0, 0).unwrap().is_empty());
    }

    #[test]
    fn read_scan_srm_v66_windows_decodes_pairs() {
        let mut bytes = 2u32.to_le_bytes().to_vec();
        bytes.extend_from_slice(&[0u8; 28]);
        bytes.extend_from_slice(&1.0f32.to_le_bytes());
        bytes.extend_from_slice(&2.0f32.to_le_bytes());
        bytes.extend_from_slice(&3.0f32.to_le_bytes());
        bytes.extend_from_slice(&4.0f32.to_le_bytes());
        let mut cursor = Cursor::new(bytes);
        let windows = read_scan_srm_v66_windows(&mut cursor, 0, 0).unwrap();
        assert_eq!(windows, vec![(1.0, 2.0), (3.0, 4.0)]);
    }

    #[test]
    fn read_flat_peaks_decodes_every_peak_of_a_multi_peak_record() {
        // The block read must decode the same (mz, abundance) pairs, in the
        // same order, that the field-at-a-time read did.
        let expected = [(500.0f32, 10.0f32), (600.5, 2.5), (700.25, 0.0)];
        let mut bytes = vec![0u8; 32]; // leading padding
        for (mz, abundance) in expected {
            bytes.extend_from_slice(&mz.to_le_bytes());
            bytes.extend_from_slice(&abundance.to_le_bytes());
        }
        // One trailing flag byte per peak: the section this reader seeks to
        // is `peak_count * 9` bytes from the end of the record.
        bytes.extend_from_slice(&[0u8; 3]);
        let cum_end = bytes.len() as u64;
        let mut cursor = Cursor::new(bytes);
        let peaks = read_flat_peaks(&mut cursor, 0, cum_end, expected.len() as u32 + 1).unwrap();
        assert_eq!(peaks.len(), expected.len());
        for (peak, (mz, abundance)) in peaks.iter().zip(expected) {
            assert_eq!((peak.mz, peak.abundance), (f64::from(mz), abundance));
        }
    }

    #[test]
    fn read_flat_peaks_rejects_a_peak_section_past_the_end_of_the_file() {
        // The scan index claims a record ending at byte 100 of a 30-byte
        // file, so the peak section starts past the end. This reader is
        // built on `BinaryReader::new`, which knows the real stream length,
        // so the allocation bound rejects the count before any byte of the
        // block is read, exactly as the field-at-a-time loop's `check_count`
        // did.
        let mut bytes = Vec::new();
        for mz in [500.0f32, 600.0] {
            bytes.extend_from_slice(&mz.to_le_bytes());
            bytes.extend_from_slice(&1.0f32.to_le_bytes());
        }
        bytes.extend_from_slice(&700.0f32.to_le_bytes()); // half of peak three
        assert_eq!(bytes.len(), 20);
        let mut cursor = Cursor::new(bytes);
        let error = read_flat_peaks(&mut cursor, 0, 100, 4).unwrap_err();
        assert!(
            matches!(error, Error::AllocationTooLarge { .. }),
            "got {error:?}"
        );
    }

    #[test]
    fn read_scan_srm_v66_decodes_every_peak_of_a_multi_peak_record() {
        let expected = [(105.0f32, 42.0f32), (205.5, 1.5), (305.25, 0.0)];
        let mut bytes = (expected.len() as u32).to_le_bytes().to_vec();
        bytes.extend_from_slice(&[0u8; 28]);
        for index in 0..expected.len() {
            // m/z window table: one (lo, hi) pair per channel.
            bytes.extend_from_slice(&(index as f32).to_le_bytes());
            bytes.extend_from_slice(&(index as f32 + 1.0).to_le_bytes());
        }
        for (index, (mz, abundance)) in expected.iter().enumerate() {
            bytes.extend_from_slice(&(index as u32).to_le_bytes()); // channel
            bytes.extend_from_slice(&mz.to_le_bytes());
            bytes.extend_from_slice(&abundance.to_le_bytes());
        }
        let mut cursor = Cursor::new(bytes);
        let peaks = read_scan_srm_v66(&mut cursor, 0, 0, 0).unwrap();
        assert_eq!(peaks.len(), expected.len());
        for (peak, (mz, abundance)) in peaks.iter().zip(expected) {
            assert_eq!((peak.mz, peak.abundance), (f64::from(mz), abundance));
        }
    }

    #[test]
    fn read_scan_srm_v66_rejects_a_record_that_ends_mid_peak() {
        // Two transitions declared, one peak record and a third of the next.
        let mut bytes = 2u32.to_le_bytes().to_vec();
        bytes.extend_from_slice(&[0u8; 28]);
        bytes.extend_from_slice(&[0u8; 16]); // two (lo, hi) window pairs
        bytes.extend_from_slice(&0u32.to_le_bytes());
        bytes.extend_from_slice(&105.0f32.to_le_bytes());
        bytes.extend_from_slice(&42.0f32.to_le_bytes());
        bytes.extend_from_slice(&1u32.to_le_bytes()); // channel of peak two only
        let mut cursor = Cursor::new(bytes);
        let error = read_scan_srm_v66(&mut cursor, 0, 0, 0).unwrap_err();
        assert!(
            matches!(error, Error::AllocationTooLarge { .. }),
            "got {error:?}"
        );
    }

    #[test]
    fn read_scan_srm_v66_windows_rejects_a_table_that_ends_mid_pair() {
        // Two pairs declared, one pair and half of the next present.
        let mut bytes = 2u32.to_le_bytes().to_vec();
        bytes.extend_from_slice(&[0u8; 28]);
        bytes.extend_from_slice(&1.0f32.to_le_bytes());
        bytes.extend_from_slice(&2.0f32.to_le_bytes());
        bytes.extend_from_slice(&3.0f32.to_le_bytes());
        let mut cursor = Cursor::new(bytes);
        let error = read_scan_srm_v66_windows(&mut cursor, 0, 0).unwrap_err();
        assert!(
            matches!(error, Error::AllocationTooLarge { .. }),
            "got {error:?}"
        );
    }

    #[test]
    fn search_v63_transition_finds_matching_record() {
        // 72-byte-per-channel table: [flag][unk][Q1][Q3_center][Q3_width][dwell][CE]...
        let mut data = vec![0u8; 8]; // leading offset so `j - 8` is in range
        data.extend_from_slice(&1.0f64.to_le_bytes()); // flag
        data.extend_from_slice(&0.0f64.to_le_bytes()); // unknown
        data.extend_from_slice(&500.0f64.to_le_bytes()); // Q1
        data.extend_from_slice(&300.0f64.to_le_bytes()); // Q3 center
        data.extend_from_slice(&1.0f64.to_le_bytes()); // Q3 width
        data.extend_from_slice(&0.02f64.to_le_bytes()); // dwell
        data.extend_from_slice(&25.0f64.to_le_bytes()); // CE
                                                        // Trailing padding: the scan only considers positions `j` with
                                                        // `j < data.len() - 32`, so there must be slack past the CE field
                                                        // for the Q3-center position (at offset 32) to be in range.
        data.extend_from_slice(&[0u8; 16]);
        let (q1, q3w, ce) = search_v63_transition(&data, 300.0).unwrap();
        assert_eq!(q1, 500.0);
        assert_eq!(q3w, 1.0);
        assert_eq!(ce, 25.0);
    }

    #[test]
    fn search_v63_transition_no_match_returns_none() {
        let data = vec![0u8; 128];
        assert!(search_v63_transition(&data, 300.0).is_none());
    }

    #[test]
    fn to_mz_intensity_capacity_uses_actual_signal_length_not_untrusted_nbins() {
        // `nbins` deliberately does not match the real chunk data below; the
        // conversion must not use it to size an allocation (see the fix in
        // `Profile::to_mz_intensity`), and the result must reflect the real
        // (small) amount of signal data actually present.
        let profile = Profile {
            first_value: 100.0,
            step: 0.01,
            peak_count: 1,
            nbins: u32::MAX,
            chunks: vec![ProfileChunk {
                first_bin: 0,
                signal: vec![1.0, 2.0, 3.0],
                fudge: None,
            }],
        };
        let result = profile.to_mz_intensity(&[]);
        assert_eq!(result.len(), 3);
    }

    #[test]
    fn to_mz_intensity_applies_fudge_in_mz_domain_after_conversion() {
        // Orbitrap-style 7-coefficient conversion: mz = A + B/f^2 + C/f^4.
        // With A = 0, C = 0 and B = 1e6, f = 1000 gives mz = 1.0 exactly.
        let coeffs = [0.0, 0.0, 0.0, 1.0e6, 0.0, 0.0, 0.0];
        let fudge = 0.25_f32;
        let profile = Profile {
            first_value: 1000.0,
            step: 0.0,
            peak_count: 1,
            nbins: 1,
            chunks: vec![ProfileChunk {
                first_bin: 0,
                signal: vec![1.0],
                fudge: Some(fudge),
            }],
        };
        let (mz, _) = profile.to_mz_intensity(&coeffs)[0];
        // Correct: 1.0 + 0.25. Frequency-domain application would instead
        // give 1e6 / 1000.25^2 = 0.99950..., a shift in the wrong direction.
        assert!((mz - 1.25).abs() < 1e-12, "got {mz}");
    }
}
