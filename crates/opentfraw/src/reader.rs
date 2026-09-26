use std::collections::HashMap;
use std::io::{Read, Seek, SeekFrom};

use crate::binary_reader::BinaryReader;
use crate::error::{Error, Result};
use crate::error_log::ErrorEntry;
use crate::generic_data::{GenericDataHeader, GenericRecord};
use crate::header::FileHeader;
use crate::metadata::{ScanParams, StatusLogEntry};
use crate::raw_file_info::RawFileInfo;
use crate::run_header::RunHeader;
use crate::scan_data::{
    read_flat_peaks, read_scan_srm_v66, search_v63_transition, Peak, ScanDataPacket,
};
use crate::scan_event::{ScanEvent, ScanEventPreamble};
use crate::scan_index::ScanIndexEntry;
use crate::seq_row::SeqRow;

/// Compute the number of scans from the run header's declared first/last
/// scan numbers.
///
/// Both inputs are untrusted u32 values read straight from the file, so
/// `last_scan - first_scan` can be as large as `u32::MAX` (e.g.
/// `first_scan == 0`, `last_scan == u32::MAX`); the `+ 1` therefore
/// saturates rather than wrapping or panicking. Any resulting oversized
/// value is caught by `BinaryReader::check_count` before it can drive an
/// allocation.
fn compute_num_scans(first_scan: u32, last_scan: u32) -> u32 {
    if last_scan >= first_scan {
        (last_scan - first_scan).saturating_add(1)
    } else {
        0
    }
}

/// Convert a 1-based scan number to a `scan_index` array offset, or `None`
/// if `scan_number` is below `first_scan_number` (a plain subtraction here
/// would underflow rather than simply being "out of range" - both
/// `scan_number` and `first_scan_number` can be arbitrary/untrusted values).
fn scan_number_to_index(scan_number: u32, first_scan_number: u32) -> Option<usize> {
    scan_number
        .checked_sub(first_scan_number)
        .map(|v| v as usize)
}

/// A parsed Thermo Fisher RAW file.
pub struct RawFileReader {
    pub header: FileHeader,
    pub seq_row: SeqRow,
    pub raw_file_info: RawFileInfo,
    pub run_header: RunHeader,
    pub scan_index: Vec<ScanIndexEntry>,
    pub scan_events: Vec<ScanEvent>,
    pub scan_parameters_header: GenericDataHeader,
    pub scan_parameters: Vec<GenericRecord>,
    pub error_log: Vec<ErrorEntry>,
    // Instrument log uses same structure
    pub inst_log_header: GenericDataHeader,
    pub inst_log: Vec<GenericRecord>,
    /// Raw file version from the header.
    pub version: u32,
    /// Number of scans.
    pub num_scans: u32,
    /// Data stream base address (for computing absolute scan offsets).
    pub data_addr: u64,
    /// True if scan data uses flat-peak format (TSQ/SRM) instead of PacketHeader.
    pub flat_peaks: bool,
    /// Detected scan-data encoding (the format used by [`Self::read_scan_peaks`]).
    pub scan_format: crate::scan_format::ScanDataFormat,
    /// Detected device family (informational).
    pub device_family: crate::device::DeviceFamily,
    /// Canonical instrument model name if one was detected in the file's
    /// metadata region (e.g. `"Orbitrap Fusion Lumos"`). `None` means only
    /// the coarse family could be inferred.
    pub instrument_model: Option<&'static str>,
    /// For SRM (flat-peak) files: maps scan_event index → Q1 precursor mass (m/z).
    ///
    /// Populated at open time by scanning the method/transition table stored
    /// in the pre-scan-data header region. Empty for non-SRM instruments.
    pub srm_q1_by_event: HashMap<u16, f64>,
    /// For SRM (flat-peak) files: maps scan_event index → Q3 isolation window pairs (lo, hi) in m/z.
    ///
    /// Populated at open time by reading the Q3 window table from the first scan record
    /// of each unique scan event class. Empty for non-SRM instruments.
    pub srm_q3_windows: HashMap<u16, Vec<(f32, f32)>>,
    /// For SRM (flat-peak) files: maps scan_event index → collision energy (eV).
    ///
    /// Populated from the v63 transition table at open time.  For v66 files the
    /// collision energy is read from per-scan parameters instead, so this map is
    /// empty for v66/TSQ Altis files.
    pub srm_ce_by_event: HashMap<u16, f64>,
}

// -- Multi-controller metadata --

/// Controller type codes as used in Thermo RAW files.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ControllerType {
    Ms,
    Analog,
    Adc,
    Pda,
    Uv,
    Other,
}

impl ControllerType {
    fn from_nsegs_ntrailer(ntrailer: u32, nsegs: u32) -> Self {
        // Heuristic: MS controller always has ntrailer > 0 (v64+) or nsegs > 0.
        // Non-MS controllers (UV, analog, PDA) have ntrailer == 0 and nsegs == 1.
        // We can't reliably distinguish between non-MS types without parsing
        // the InstID/method block, so we fall back to Other for those.
        if ntrailer > 0 || nsegs > 1 {
            Self::Ms
        } else {
            Self::Other
        }
    }
}

/// Minimal metadata about one controller in a multi-controller RAW file.
#[derive(Debug, Clone)]
pub struct ControllerInfo {
    /// Zero-based controller index (position in `run_header_addrs`).
    pub index: usize,
    /// File offset to this controller's RunHeader.
    pub run_header_addr: u64,
    /// Whether this controller is the primary MS controller.
    pub is_ms_controller: bool,
    /// Inferred controller type.
    pub controller_type: ControllerType,
    /// First scan number.
    pub first_scan: u32,
    /// Last scan number.
    pub last_scan: u32,
    /// Acquisition start time (minutes).
    pub start_time: f64,
    /// Acquisition end time (minutes).
    pub end_time: f64,
}

impl RawFileReader {
    /// Open and parse a RAW file from a reader.
    pub fn open<R: Read + Seek>(source: R) -> Result<Self> {
        Self::open_inner(source)
            .map_err(|error| error.with_context(crate::ErrorContext::new(crate::Operation::Open)))
    }

    fn open_inner<R: Read + Seek>(source: R) -> Result<Self> {
        let mut r = BinaryReader::new(source);

        // 1. FileHeader
        let header = FileHeader::read(&mut r)?;
        let version = header.version;

        // 2. SeqRow
        let seq_row = SeqRow::read(&mut r, version)?;

        // 3. ASInfo (read and discard preamble + string)
        let _as_preamble = r.read_bytes(24)?; // ASInfoPreamble: 24 bytes
        let _as_text = r.read_pascal_string()?;

        // 4. RawFileInfo
        let raw_file_info = RawFileInfo::read(&mut r, version)?;

        // 5. Extract addresses
        let data_addr = raw_file_info.preamble.data_addr;

        // 6. Select the MS controller RunHeader.
        // Multi-controller files (e.g. UV + MS) have one RunHeader per controller.
        // The MS controller has ntrailer > 0 (v64+) or first_scan <= last_scan with
        // nsegs > 0 (v63 and earlier). We iterate all addresses and pick the best.
        let run_header = {
            let addrs = &raw_file_info.preamble.run_header_addrs;
            let mut chosen = None;
            for &addr in addrs {
                if addr == 0 {
                    continue;
                }
                r.seek_to(addr)?;
                let rh = RunHeader::read(&mut r, version)?;
                // Heuristic for identifying the MS controller:
                // 1. For v64+: ntrailer > 0 (scan events present) - catches most instruments.
                // 2. For all versions: RunHeader.data_addr == preamble.data_addr - the MS
                //    controller's scan data begins at the same address the preamble declares.
                //    This catches TSQ/triple-quad instruments where ntrailer=0 (no scan events).
                // 3. Pre-v64 fallback: valid scan range with nsegs > 0.
                let is_ms = if version >= 64 {
                    rh.ntrailer > 0 || rh.data_addr == data_addr
                } else {
                    rh.sample_info.last_scan_number >= rh.sample_info.first_scan_number
                        && rh.nsegs > 0
                };
                if is_ms {
                    chosen = Some(rh);
                    break;
                }
            }
            // Fall back to first address if no MS controller found
            match chosen {
                Some(rh) => rh,
                None => {
                    r.seek_to(addrs[0])?;
                    RunHeader::read(&mut r, version)?
                }
            }
        };

        let first_scan = run_header.sample_info.first_scan_number;
        let last_scan = run_header.sample_info.last_scan_number;
        let num_scans = compute_num_scans(first_scan, last_scan);

        // 7. Scan index
        r.seek_to(run_header.scan_index_addr)?;
        r.check_count(
            num_scans as u64,
            ScanIndexEntry::min_size_for_version(version),
        )?;
        let mut scan_index = Vec::with_capacity(num_scans as usize);
        for _ in 0..num_scans {
            scan_index.push(ScanIndexEntry::read(&mut r, version)?);
        }

        // 8. Scan event trailer
        r.seek_to(run_header.scan_trailer_addr)?;
        let n_events = if version >= 64 {
            // v64+: first u32 is a preamble (not count); use ntrailer from RunHeader
            let _preamble = r.read_u32()?;
            run_header.ntrailer
        } else {
            r.read_u32()?
        };
        let scan_events = if let Some(events) =
            crate::qualified_events::read(&mut r, version, n_events, run_header.scan_params_addr)?
        {
            events
        } else {
            // For v66, compute per-event body sizes from the stream's address range.
            // The scan event stream spans [scan_trailer_addr+4 .. scan_params_addr).
            // Each event = preamble (136 bytes) + body.
            //
            // Simple instruments (Q Exactive, Exploris): all events are identical in
            // size so stream_bytes divides evenly by n_events.
            //
            // Tribrid instruments (Eclipse, Fusion Lumos): primary (MS1) scans and
            // dependent (MS2+) scans have different body layouts:
            //   Primary event:   232 bytes total (preamble 136 + body 96)
            //   Dependent event: 344 bytes total (preamble 136 + body 208)
            // Confirmed empirically across Orbitrap Eclipse (EThcD) and Fusion Lumos
            // (DIA, MS3) files.
            let preamble_size = ScanEventPreamble::size_for_version(version);
            let (v66_body_primary, v66_body_dependent): (usize, usize) = if version >= 66
                && n_events > 0
            {
                let stream_bytes = run_header
                    .scan_params_addr
                    .saturating_sub(run_header.scan_trailer_addr)
                    .saturating_sub(4);
                let remainder = stream_bytes % n_events as u64;
                if remainder == 0 {
                    // Uniform event size (Q Exactive, Exploris, etc.)
                    let body = (stream_bytes / n_events as u64) as usize;
                    let body = body.saturating_sub(preamble_size);
                    (body, body)
                } else {
                    // Variable-length events: tribrid Orbitrap instruments.
                    // Known sizes: primary=232, dependent=344 (body 96 and 208).
                    const PRIMARY_EVENT: u64 = 232;
                    const DEPENDENT_EVENT: u64 = 344;
                    let gap = DEPENDENT_EVENT - PRIMARY_EVENT;
                    let n = n_events as u64;
                    // n_primary * PRIMARY_EVENT + n_dependent * DEPENDENT_EVENT = stream_bytes
                    // n_primary + n_dependent = n
                    // => n_primary = (n * DEPENDENT_EVENT - stream_bytes) / gap
                    let n_primary_numerator = n
                        .saturating_mul(DEPENDENT_EVENT)
                        .saturating_sub(stream_bytes);
                    if n_primary_numerator % gap == 0 {
                        let n_primary = n_primary_numerator / gap;
                        let n_dependent = n.saturating_sub(n_primary);
                        let total_check = n_primary * PRIMARY_EVENT + n_dependent * DEPENDENT_EVENT;
                        if total_check == stream_bytes {
                            // Verified: use the tribrid sizes.
                            (
                                (PRIMARY_EVENT as usize).saturating_sub(preamble_size),
                                (DEPENDENT_EVENT as usize).saturating_sub(preamble_size),
                            )
                        } else {
                            return Err(Error::UnsupportedLayout("variable event layout"));
                        }
                    } else {
                        return Err(Error::UnsupportedLayout("variable event layout"));
                    }
                }
            } else {
                (0, 0)
            };
            // Each event consumes at least its preamble on disk (the body may add
            // more, but the preamble alone is a safe lower bound for the count check).
            r.check_count(n_events as u64, preamble_size as u64)?;
            let mut scan_events = Vec::with_capacity(n_events as usize);
            for _ in 0..n_events {
                scan_events.push(ScanEvent::read(
                    &mut r,
                    version,
                    v66_body_primary,
                    v66_body_dependent,
                )?);
            }

            scan_events
        };

        // 9. Error log
        let n_errors = run_header.sample_info.error_log_length;
        let error_log = if n_errors > 0 {
            r.seek_to(run_header.error_log_addr)?;
            if version >= 64 {
                let _preamble = r.read_u32()?;
            }
            // Minimum on-disk size of one entry: f32 time (4) + u32 pascal
            // string char count (4), i.e. an empty message.
            r.check_count(n_errors as u64, 8)?;
            let mut log = Vec::with_capacity(n_errors as usize);
            for _ in 0..n_errors {
                log.push(ErrorEntry::read(&mut r)?);
            }
            log
        } else {
            // Ensure reader is positioned at error_log_addr even when empty.
            r.seek_to(run_header.error_log_addr)?;
            Vec::new()
        };
        // The GDH for scan parameters immediately follows the error-log entries.
        // Do NOT seek back to error_log_addr - doing so would cause find_forward
        // to scan over the scan_index (which may sit between error_log and
        // scan_trailer in some file layouts), creating a CPU-spinning O(n) search
        // through megabytes of binary scan data.
        let after_error_log = r.position();

        // 10. Scan parameters (trailer extra) - GenericData format in v64+.
        //     The schema (GDH) is written just after the error-log entries;
        //     the records are written at `scan_params_addr` (tail of file)
        //     with NO stream preamble - records begin directly at
        //     scan_params_addr. Any bytes after the last record are trailing
        //     padding and can be ignored.
        let (scan_parameters_header, scan_parameters) = if version >= 64 {
            // Search from after the error log entries up to scan_trailer.
            // This skips any scan_index data that may sit in between.
            let scan_distance = run_header.scan_trailer_addr.saturating_sub(after_error_log);
            // Estimate per-record size from the tail of the file using integer
            // division. Any remainder bytes are trailing data, not a preamble.
            let file_size = r.length()?;
            let tail = file_size.saturating_sub(run_header.scan_params_addr);
            let expected_record_size = if num_scans > 0 && tail > 0 {
                let per_scan = tail / num_scans as u64;
                if per_scan >= 4 {
                    Some(per_scan as usize)
                } else {
                    None
                }
            } else {
                None
            };
            match GenericDataHeader::find_forward(&mut r, scan_distance, expected_record_size)? {
                Some(hdr) => {
                    // Records start directly at scan_params_addr - no stream preamble.
                    r.seek_to(run_header.scan_params_addr)?;
                    r.check_count(num_scans as u64, hdr.fixed_record_size().max(1) as u64)?;
                    let mut params = Vec::with_capacity(num_scans as usize);
                    for _ in 0..num_scans {
                        params.push(GenericRecord::read(&mut r, &hdr)?);
                    }
                    (hdr, params)
                }
                None => (GenericDataHeader { fields: Vec::new() }, Vec::new()),
            }
        } else {
            (GenericDataHeader { fields: Vec::new() }, Vec::new())
        };

        // 11. Instrument log - GenericData format in v64+
        let (inst_log_header, inst_log) = if version >= 64 {
            r.seek_to(run_header.inst_log_addr)?;
            match GenericDataHeader::try_read(&mut r)? {
                Some(hdr) => {
                    let n_inst = run_header.sample_info.inst_log_length;
                    r.check_count(n_inst as u64, hdr.fixed_record_size().max(1) as u64)?;
                    let mut log = Vec::with_capacity(n_inst as usize);
                    for _ in 0..n_inst {
                        log.push(GenericRecord::read(&mut r, &hdr)?);
                    }
                    (hdr, log)
                }
                None => (GenericDataHeader { fields: Vec::new() }, Vec::new()),
            }
        } else {
            (GenericDataHeader { fields: Vec::new() }, Vec::new())
        };

        // Detect flat-peak (TSQ/SRM) format.
        // Reliable indicator: ntrailer == 0 means no scan event trailer was written, which
        // is the case for all TSQ/triple-quad SRM instruments.
        // Fallback: first scan data_size < 100 (catches edge cases with tiny SRM windows).
        // In the flat format, data_size is the number of MRM peaks, not bytes.
        let flat_peaks = run_header.ntrailer == 0
            || scan_index
                .first()
                .map(|e| e.data_size < 100)
                .unwrap_or(false);

        // Classify scan format and device family.
        let scan_format = crate::scan_format::ScanDataFormat::detect(version, flat_peaks);
        let first_analyzer = scan_events.first().and_then(|e| e.preamble.analyzer());

        // For SRM (flat-peak) files, read the entire pre-scan-data region so that
        // we can extract Q1 values from the method/transition table stored there.
        // For other instruments, read only 64 KB for instrument model detection.
        let scan_window_cap = if flat_peaks { data_addr } else { 64 * 1024u64 };
        let window_len = scan_window_cap.min(data_addr);
        let metadata_window = if window_len > 0 {
            r.seek_to(0)?;
            r.read_bytes(window_len as usize).unwrap_or_default()
        } else {
            Vec::new()
        };
        // All BinaryReader operations are complete; reclaim the underlying source so
        // it can be used for on-demand reads (e.g. Q3 window table from scan records).
        let mut source = r.into_inner();
        let detected = crate::device::DeviceFamily::detect_instrument(
            &metadata_window,
            &header.audit_start.tag2,
            &seq_row.inst_method,
            first_analyzer,
        );
        let device_family = detected.family;
        let instrument_model = detected.model;

        // For SRM files: extract Q1 masses, Q3 window pairs, and (for v63) collision energies
        // from the pre-scan-data header region and/or the scan data records.
        //
        // v66 (TSQ Quantiva / TSQ Altis, FlatV66):
        //   Transition table layout: [Q1: f64][Q3_lo: f64][Q3_hi: f64] per channel.
        //   Anchor: scan_index.high_mz equals the Q3_hi of the highest-Q3 channel for each
        //   event class.  Q3 window pairs come from the per-scan record header.
        //
        // v63 (TSQ Quantum / TSQ Vantage, FlatV63):
        //   Transition table layout: 72-byte records; Q1 at [+16], Q3_center at [+24],
        //   Q3_width at [+32], CE at [+48].  scan_index.low_mz/high_mz hold the instrument
        //   scan range (not per-transition values), so the high_mz anchor does not apply.
        //   Q3 centers come from the first scan's peak list; Q3 windows are computed as
        //   Q3_center ± Q3_width/2.
        let (srm_q1_by_event, srm_q3_windows, srm_ce_by_event) = {
            use crate::scan_format::ScanDataFormat;
            match (flat_peaks, scan_format) {
                (true, ScanDataFormat::FlatV66) if metadata_window.len() >= 24 => {
                    // --- v66 Q1 extraction: anchor on scan_index.high_mz ---
                    let mut event_q3_hi: HashMap<u16, f64> = HashMap::new();
                    for entry in &scan_index {
                        if entry.high_mz > 50.0 && entry.high_mz < 2000.0 {
                            event_q3_hi.entry(entry.scan_event).or_insert(entry.high_mz);
                        }
                    }
                    let data = &metadata_window;
                    let mut q1_map: HashMap<u16, f64> = HashMap::new();
                    'outer_v66: for (&event, &q3_hi_target) in &event_q3_hi {
                        let end = data.len().saturating_sub(8);
                        for i in 16..end {
                            let hi = crate::bytes::read_f64_le(data, i)?;
                            if (hi - q3_hi_target).abs() < 0.002 {
                                let lo = crate::bytes::read_f64_le(data, i - 8)?;
                                if hi > lo && (hi - lo) < 0.1 {
                                    let q1 = crate::bytes::read_f64_le(data, i - 16)?;
                                    if q1 > 50.0 && q1 < 3000.0 {
                                        q1_map.insert(event, q1);
                                        continue 'outer_v66;
                                    }
                                }
                            }
                        }
                    }
                    // --- v66 Q3 window extraction: read per-scan record header ---
                    let mut seen: HashMap<u16, bool> = HashMap::new();
                    let mut q3_map: HashMap<u16, Vec<(f32, f32)>> = HashMap::new();
                    for entry in &scan_index {
                        if seen.contains_key(&entry.scan_event) {
                            continue;
                        }
                        seen.insert(entry.scan_event, true);
                        if let Ok(windows) = crate::scan_data::read_scan_srm_v66_windows(
                            &mut source,
                            data_addr,
                            entry.offset,
                        ) {
                            if !windows.is_empty() {
                                q3_map.insert(entry.scan_event, windows);
                            }
                        }
                    }
                    (q1_map, q3_map, HashMap::new())
                }
                (true, ScanDataFormat::FlatV63) => {
                    // --- v63 Q1 + Q3 window + CE extraction ---
                    // Read peaks from the first scan of each event class; each peak's mz
                    // is the Q3 center for that channel.  Search the pre-data region for
                    // the Q3_center value to find Q1, Q3_width, and CE from the transition
                    // table.  Q3 windows are computed as (Q3_center - width/2, Q3_center + width/2).
                    let mut seen: HashMap<u16, bool> = HashMap::new();
                    let mut q1_map: HashMap<u16, f64> = HashMap::new();
                    let mut q3_map: HashMap<u16, Vec<(f32, f32)>> = HashMap::new();
                    let mut ce_map: HashMap<u16, f64> = HashMap::new();
                    let data = &metadata_window;
                    for entry in &scan_index {
                        let ev = entry.scan_event;
                        if seen.contains_key(&ev) {
                            continue;
                        }
                        seen.insert(ev, true);
                        let peaks = match read_flat_peaks(
                            &mut source,
                            data_addr,
                            entry.offset,
                            entry.data_size,
                        ) {
                            Ok(p) if !p.is_empty() => p,
                            _ => continue,
                        };
                        // Use the first peak's mz as Q3_center anchor.
                        if let Some((q1, q3w, ce)) = search_v63_transition(data, peaks[0].mz) {
                            q1_map.insert(ev, q1);
                            ce_map.insert(ev, ce);
                            let half = (q3w / 2.0) as f32;
                            let windows: Vec<(f32, f32)> = peaks
                                .iter()
                                .map(|p| (p.mz as f32 - half, p.mz as f32 + half))
                                .collect();
                            q3_map.insert(ev, windows);
                        }
                    }
                    (q1_map, q3_map, ce_map)
                }
                _ => (HashMap::new(), HashMap::new(), HashMap::new()),
            }
        };

        Ok(Self {
            header,
            seq_row,
            raw_file_info,
            run_header,
            scan_index,
            scan_events,
            scan_parameters_header,
            scan_parameters,
            error_log,
            inst_log_header,
            inst_log,
            version,
            num_scans,
            data_addr,
            flat_peaks,
            scan_format,
            device_family,
            instrument_model,
            srm_q1_by_event,
            srm_q3_windows,
            srm_ce_by_event,
        })
    }

    /// Open a RAW file from a path.
    pub fn open_path(path: impl AsRef<std::path::Path>) -> Result<Self> {
        let path = path.as_ref();
        let result = (|| {
            let file = std::fs::File::open(path)?;
            Self::open(std::io::BufReader::new(file))
        })();
        result.map_err(|error: Error| {
            let mut context = crate::ErrorContext::new(crate::Operation::Open);
            context.path = Some(path.to_path_buf());
            error.with_context(context)
        })
    }

    /// Enumerate all controllers in this RAW file.
    ///
    /// Multi-detector acquisition systems write one [`RunHeader`] per
    /// controller (MS, UV, PDA, Analog). This method parses all controller
    /// headers and returns a `Vec<ControllerInfo>` with basic metadata for
    /// each. The primary MS controller can be identified via
    /// [`ControllerInfo::is_ms_controller`].
    ///
    /// For single-controller files (the common case), this returns a
    /// one-element vec with the MS controller.
    pub fn controllers<R: Read + Seek>(&self, source: &mut R) -> Result<Vec<ControllerInfo>> {
        let mut r = BinaryReader::new(source);
        let addrs = &self.raw_file_info.preamble.run_header_addrs;
        let mut infos = Vec::with_capacity(addrs.len());
        for (i, &addr) in addrs.iter().enumerate() {
            if addr == 0 {
                continue;
            }
            r.seek_to(addr)?;
            let rh = RunHeader::read(&mut r, self.version)?;
            let is_ms = if self.version >= 64 {
                rh.ntrailer > 0 || rh.data_addr == self.data_addr
            } else {
                rh.nsegs > 0
            };
            let ct = if is_ms {
                ControllerType::Ms
            } else {
                ControllerType::from_nsegs_ntrailer(rh.ntrailer, rh.nsegs)
            };
            infos.push(ControllerInfo {
                index: i,
                run_header_addr: addr,
                is_ms_controller: is_ms,
                controller_type: ct,
                first_scan: rh.sample_info.first_scan_number,
                last_scan: rh.sample_info.last_scan_number,
                start_time: rh.sample_info.start_time,
                end_time: rh.sample_info.end_time,
            });
        }
        Ok(infos)
    }

    /// Look up the [`ScanIndexEntry`] for a 1-based scan number, bounds- and
    /// underflow-checked.
    ///
    /// `scan_number` is caller-provided (and, via the public API, ultimately
    /// file/user controlled) and `first_scan_number` is itself an untrusted
    /// value read from the file, so a plain `scan_number - first_scan_number`
    /// subtraction can underflow when `scan_number` is smaller - not just an
    /// "index too large" situation. Both cases are simply out of range.
    fn scan_index_entry(&self, scan_number: u32) -> Result<&ScanIndexEntry> {
        let idx = scan_number_to_index(scan_number, self.run_header.sample_info.first_scan_number)
            .ok_or(Error::InvalidScanNumber(scan_number))?;
        self.scan_index
            .get(idx)
            .ok_or(Error::InvalidScanNumber(scan_number))
    }

    pub(crate) fn scan_context(
        &self,
        scan_number: u32,
        operation: crate::Operation,
    ) -> crate::ErrorContext {
        let scan_index = scan_number.checked_sub(self.run_header.sample_info.first_scan_number);
        let byte_offset = self.scan_index_entry(scan_number).ok().and_then(|entry| {
            let offset = self.data_addr.checked_add(entry.offset)?;
            if self.scan_format == crate::ScanDataFormat::FlatV63 {
                offset.checked_sub(u64::from(entry.data_size))
            } else {
                Some(offset)
            }
        });
        crate::ErrorContext {
            operation,
            path: None,
            scan_index,
            scan_number: Some(scan_number),
            byte_offset,
        }
    }

    /// Read a single scan data packet (PacketHeader format).
    pub fn read_scan<R: Read + Seek>(
        &self,
        source: &mut R,
        scan_number: u32,
    ) -> Result<ScanDataPacket> {
        self.read_scan_inner(source, scan_number).map_err(|error| {
            error.with_context(self.scan_context(scan_number, crate::Operation::ReadScan))
        })
    }

    fn read_scan_inner<R: Read + Seek>(
        &self,
        source: &mut R,
        scan_number: u32,
    ) -> Result<ScanDataPacket> {
        let entry = self.scan_index_entry(scan_number)?;
        // Both are untrusted u64 values read from the file (a
        // RawFileInfo address and a ScanIndexEntry offset); saturate
        // instead of panicking on overflow for a corrupt combination.
        let abs_offset = self.data_addr.saturating_add(entry.offset);
        source.seek(SeekFrom::Start(abs_offset))?;
        let mut r = BinaryReader::new(source);
        ScanDataPacket::read(&mut r)
    }

    /// Read a single scan packet's centroid peaks and FT label data
    /// (resolution / noise / baseline), skipping the profile signal for speed.
    ///
    /// Only valid for PacketHeader-format files; TSQ/SRM scans carry no FT
    /// label data.
    pub fn read_scan_labels<R: Read + Seek>(
        &self,
        source: &mut R,
        scan_number: u32,
    ) -> Result<ScanDataPacket> {
        self.read_scan_labels_inner(source, scan_number)
            .map_err(|error| {
                error.with_context(self.scan_context(scan_number, crate::Operation::ReadScan))
            })
    }

    fn read_scan_labels_inner<R: Read + Seek>(
        &self,
        source: &mut R,
        scan_number: u32,
    ) -> Result<ScanDataPacket> {
        use crate::scan_format::ScanDataFormat;
        if self.scan_format != ScanDataFormat::PacketHeader {
            return Err(Error::UnsupportedOperation(
                "centroid_labels / read_scan_labels requires a PacketHeader file (Orbitrap/ion-trap); TSQ/SRM files carry no FT label data",
            ));
        }
        let entry = self.scan_index_entry(scan_number)?;
        // Both are untrusted u64 values read from the file (a
        // RawFileInfo address and a ScanIndexEntry offset); saturate
        // instead of panicking on overflow for a corrupt combination.
        let abs_offset = self.data_addr.saturating_add(entry.offset);
        source.seek(SeekFrom::Start(abs_offset))?;
        let mut r = BinaryReader::new(source);
        ScanDataPacket::read_skip_profile(&mut r)
    }

    /// Read a single scan as flat peaks (TSQ/SRM format).
    ///
    /// In this format, `entry.offset` is the cumulative end byte offset within
    /// the data stream. Peaks are (f32, f32) pairs at the end of each record.
    pub fn read_scan_flat<R: Read + Seek>(
        &self,
        source: &mut R,
        scan_number: u32,
    ) -> Result<Vec<Peak>> {
        self.read_scan_flat_inner(source, scan_number)
            .map_err(|error| {
                error.with_context(self.scan_context(scan_number, crate::Operation::ReadScan))
            })
    }

    fn read_scan_flat_inner<R: Read + Seek>(
        &self,
        source: &mut R,
        scan_number: u32,
    ) -> Result<Vec<Peak>> {
        let entry = self.scan_index_entry(scan_number)?;
        read_flat_peaks(source, self.data_addr, entry.offset, entry.data_size)
    }

    /// Read a single scan in v66 SRM format (TSQ Quantiva / TSQ Altis).
    ///
    /// `entry.offset` is the START byte offset within the data stream.
    /// The record is fixed-size (`entry.data_size` bytes) and contains:
    ///   n_peaks (u32), header, m/z window table, then peak triplets.
    pub fn read_scan_srm_v66<R: Read + Seek>(
        &self,
        source: &mut R,
        scan_number: u32,
    ) -> Result<Vec<Peak>> {
        self.read_scan_srm_v66_inner(source, scan_number)
            .map_err(|error| {
                error.with_context(self.scan_context(scan_number, crate::Operation::ReadScan))
            })
    }

    fn read_scan_srm_v66_inner<R: Read + Seek>(
        &self,
        source: &mut R,
        scan_number: u32,
    ) -> Result<Vec<Peak>> {
        let entry = self.scan_index_entry(scan_number)?;
        read_scan_srm_v66(source, self.data_addr, entry.offset, entry.data_size)
    }

    /// Read a single scan's peaks using whichever decoder matches this file's
    /// scan-data format.
    ///
    /// This low-level entry point dispatches on
    /// [`Self::scan_format`] so callers do not have to know whether a file is
    /// a TSQ SRM run (flat peaks) or an Orbitrap/ion-trap acquisition
    /// (PacketHeader records). Prefer [`crate::try_extract_spectrum_arrays`]
    /// for validated scientific arrays and contextual metadata.
    ///
    /// The returned `Vec<Peak>` contains centroided peaks regardless of the
    /// underlying format. For PacketHeader files that also contain a profile
    /// signal, use [`Self::read_scan`] to access both.
    pub fn read_scan_peaks<R: Read + Seek>(
        &self,
        source: &mut R,
        scan_number: u32,
    ) -> Result<Vec<Peak>> {
        use crate::scan_format::ScanDataFormat;
        match self.scan_format {
            ScanDataFormat::PacketHeader => {
                let pkt = self.read_scan(source, scan_number)?;
                Ok(pkt.peaks)
            }
            ScanDataFormat::FlatV63 => self.read_scan_flat(source, scan_number),
            ScanDataFormat::FlatV66 => self.read_scan_srm_v66(source, scan_number),
        }
    }

    /// Read centroided peaks only, skipping profile data.
    ///
    /// For PacketHeader files (Orbitrap / ion-trap), this skips the large
    /// profile-data section, making it 2-10× faster than
    /// [`Self::read_scan_peaks`] when only centroided m/z and intensity values
    /// are needed (e.g. mzML export, peak area queries).
    ///
    /// For TSQ/SRM files this is identical to [`Self::read_scan_peaks`].
    pub fn read_peaks_only<R: Read + Seek>(
        &self,
        source: &mut R,
        scan_number: u32,
    ) -> Result<Vec<Peak>> {
        self.read_peaks_only_inner(source, scan_number)
            .map_err(|error| {
                error.with_context(self.scan_context(scan_number, crate::Operation::ReadScan))
            })
    }

    fn read_peaks_only_inner<R: Read + Seek>(
        &self,
        source: &mut R,
        scan_number: u32,
    ) -> Result<Vec<Peak>> {
        use crate::scan_format::ScanDataFormat;
        match self.scan_format {
            ScanDataFormat::PacketHeader => {
                let entry = self.scan_index_entry(scan_number)?;
                // Both are untrusted u64 values read from the file (a
                // RawFileInfo address and a ScanIndexEntry offset); saturate
                // instead of panicking on overflow for a corrupt combination.
                let abs_offset = self.data_addr.saturating_add(entry.offset);
                source.seek(SeekFrom::Start(abs_offset))?;
                let mut r = BinaryReader::new(source);
                ScanDataPacket::read_peaks_only(&mut r)
            }
            ScanDataFormat::FlatV63 => self.read_scan_flat(source, scan_number),
            ScanDataFormat::FlatV66 => self.read_scan_srm_v66(source, scan_number),
        }
    }

    /// Return the scan-parameter record for a given 1-based scan number.
    ///
    /// Returns `None` if the file has no scan-parameter stream or if
    /// `scan_number` is outside the valid scan range.
    pub fn scan_parameters(&self, scan_number: u32) -> Option<&GenericRecord> {
        let first = self.run_header.sample_info.first_scan_number;
        let idx = scan_number.checked_sub(first)? as usize;
        self.scan_parameters.get(idx)
    }

    /// Return a typed view of the scan-parameter record for a given scan.
    ///
    /// This wraps [`Self::scan_parameters`] in a [`ScanParams`] accessor that
    /// provides named, type-safe fields and handles label-name variations
    /// across instrument families.
    pub fn scan_params(&self, scan_number: u32) -> Option<ScanParams<'_>> {
        self.scan_parameters(scan_number).map(ScanParams)
    }

    /// Return the raw instrument-log record for a given scan number, or
    /// `None` if the scan is out of range or no instrument log was found.
    ///
    /// The instrument log contains per-scan instrument-state values:
    /// temperatures, voltages, pressures, ion counts, etc.
    pub fn inst_log_record(&self, scan_number: u32) -> Option<&GenericRecord> {
        let first = self.run_header.sample_info.first_scan_number;
        let idx = scan_number.checked_sub(first)? as usize;
        self.inst_log.get(idx)
    }

    /// Return a typed [`StatusLogEntry`] view for the given scan number.
    ///
    /// This wraps [`Self::inst_log_record`] and provides named, type-safe
    /// accessors for common instrument-status fields.
    pub fn status_log_entry(&self, scan_number: u32) -> Option<StatusLogEntry<'_>> {
        self.inst_log_record(scan_number).map(StatusLogEntry)
    }

    /// Return a reconstructed Thermo scan filter subset for a given scan
    /// (1-based scan number), or `None` if the scan is out of range.
    ///
    /// Example output: `"FTMS + p NSI Full ms [350.0000-1500.0000]"`.
    ///
    /// Activation clauses are included only for qualified v66 reaction tags.
    /// See [`crate::scan_filter`] for grammar details.
    pub fn scan_filter(&self, scan_number: u32) -> Option<String> {
        let first = self.run_header.sample_info.first_scan_number;
        let idx = scan_number.checked_sub(first)? as usize;
        let entry = self.scan_index.get(idx)?;

        // SRM files have no scan events; build the filter string from
        // the pre-loaded Q1 and Q3 window maps.
        if self.flat_peaks {
            let q1 = self.srm_q1_by_event.get(&entry.scan_event).copied()?;
            let windows = self.srm_q3_windows.get(&entry.scan_event)?;
            // v63 (TSQ Quantum/Vantage): NSI ionization, @cid{CE:.2} after Q1.
            // v66 (TSQ Quantiva/Altis): ESI ionization, no CE in filter.
            use crate::scan_format::ScanDataFormat;
            let ionization = match self.scan_format {
                ScanDataFormat::FlatV63 => "NSI",
                _ => "ESI",
            };
            let mut s = format!("c {} SRM ms2 {:.3}", ionization, q1);
            if !windows.is_empty() {
                s.push(' ');
                s.push('[');
                for (i, (lo, hi)) in windows.iter().enumerate() {
                    if i > 0 {
                        s.push_str(", ");
                    }
                    s.push_str(&format!("{:.3}-{:.3}", lo, hi));
                }
                s.push(']');
            }
            return Some(s);
        }

        let event = self.scan_events.get(idx)?;
        // Precursor m/z and activation energy come from the per-scan params
        // table (not the event body) for v66+. Fall back silently if missing.
        let params = self.scan_params(scan_number);
        let precursor = params.as_ref().and_then(|p| p.monoisotopic_mz());
        let energy = params.as_ref().and_then(|p| p.activation_energy());
        let supplemental = params
            .as_ref()
            .and_then(|p| p.supplemental_activation_energy());
        Some(crate::scan_filter::build_filter_with_activation(
            event,
            entry,
            precursor,
            energy,
            supplemental,
            event.qualified_activation(self.version),
        ))
    }

    /// Return all scan retention times (minutes) in scan order (1-based scan numbers).
    ///
    /// This is equivalent to collecting `scan_index[i].start_time` for every scan.
    /// The returned `Vec` is indexed by `scan_number - first_scan_number`.
    pub fn retention_times(&self) -> Vec<f64> {
        self.scan_index.iter().map(|e| e.start_time).collect()
    }

    /// Return a per-scan chromatogram as `(retention_time_min, tic)` pairs.
    pub fn tic_chromatogram(&self) -> Vec<(f64, f64)> {
        self.scan_index
            .iter()
            .map(|e| (e.start_time, e.total_current))
            .collect()
    }

    /// Return a per-scan base-peak chromatogram as `(retention_time_min, bpi, base_mz)` triples.
    pub fn bpc_chromatogram(&self) -> Vec<(f64, f64, f64)> {
        self.scan_index
            .iter()
            .map(|e| (e.start_time, e.base_intensity, e.base_mz))
            .collect()
    }

    /// Return the instrument method file path or name as stored in the
    /// sequence row. This is the name of the method used during acquisition
    /// (e.g. `"Standard_HCD.meth"`), not the embedded method text.
    ///
    /// See also [`Self::instrument_method_text`] for extracting the embedded
    /// XML/text method body from the file.
    pub fn instrument_method_name(&self) -> &str {
        &self.seq_row.inst_method
    }

    /// Attempt to extract the embedded instrument method text from the RAW file.
    ///
    /// Thermo RAW files embed the acquisition method as a UTF-16LE text or
    /// XML blob in the metadata region. This method scans the bytes between
    /// the start of the file and the scan data for the longest contiguous
    /// block of valid UTF-16LE text (at least 256 characters long) and returns
    /// it as a `String`.
    ///
    /// Returns `None` if no suitable text block is found or if the method was
    /// not embedded (`method_file_present == false`).
    ///
    /// Note: This is a best-effort extraction. The result is the raw text
    /// content; callers may wish to trim or parse it further.
    pub fn instrument_method_text<R: Read + Seek>(&self, source: &mut R) -> Option<String> {
        if !self.raw_file_info.preamble.method_file_present {
            return None;
        }
        // Read metadata region: from byte 0 up to (but not including) scan data.
        // Cap at 512 KB to avoid reading very large files entirely.
        const MAX_WINDOW: u64 = 512 * 1024;
        let window_len = MAX_WINDOW.min(self.data_addr) as usize;
        if window_len < 4 {
            return None;
        }
        source.seek(std::io::SeekFrom::Start(0)).ok()?;
        let mut buf = vec![0u8; window_len];
        source.read_exact(&mut buf).ok()?;

        // Scan for the longest valid UTF-16LE text block (min 256 chars = 512 bytes).
        // Strategy: find aligned 2-byte sequences where every pair decodes to a
        // printable/whitespace Unicode scalar (U+0020..U+FFFD).
        extract_utf16le_text(&buf, 256)
    }
}

/// Scan `buf` for the longest contiguous UTF-16LE text block of at least
/// `min_chars` characters and return it as a String. Returns `None` if no
/// such block exists.
fn extract_utf16le_text(buf: &[u8], min_chars: usize) -> Option<String> {
    if buf.len() < 2 {
        return None;
    }
    let mut best: Option<String> = None;
    let mut best_len = 0usize;

    // Try each even alignment (0 or 1 byte offset from start).
    for alignment in 0..2usize {
        let start = alignment;
        let usable = buf.len().saturating_sub(start);
        let n_units = usable / 2;
        if n_units < min_chars {
            continue;
        }

        let mut run_start = 0usize;
        let mut run_chars: Vec<u16> = Vec::with_capacity(min_chars);

        let flush = |run_chars: &Vec<u16>,
                     run_start: usize,
                     best: &mut Option<String>,
                     best_len: &mut usize| {
            if run_chars.len() >= min_chars {
                if let Ok(s) = String::from_utf16(run_chars) {
                    let _ = run_start; // suppress unused warning
                    if run_chars.len() > *best_len {
                        *best_len = run_chars.len();
                        *best = Some(s);
                    }
                }
            }
        };

        for i in 0..n_units {
            let off = start + i * 2;
            let u = u16::from_le_bytes([buf[off], buf[off + 1]]);
            let is_ok = matches!(u, 0x0009 | 0x000A | 0x000D | 0x0020..=0xFFFD);
            if is_ok {
                run_chars.push(u);
            } else {
                flush(&run_chars, run_start, &mut best, &mut best_len);
                run_start = i + 1;
                run_chars.clear();
            }
        }
        flush(&run_chars, run_start, &mut best, &mut best_len);
    }
    best
}

#[cfg(test)]
mod tests {
    use super::*;

    // Regression tests for a fuzzer-found crash: with `first_scan = 0` and
    // `last_scan = u32::MAX`, the old `last_scan - first_scan + 1`
    // computation overflowed u32 and panicked under debug assertions.
    #[test]
    fn compute_num_scans_normal_range() {
        assert_eq!(compute_num_scans(1, 5), 5);
        assert_eq!(compute_num_scans(1, 1), 1);
    }

    #[test]
    fn compute_num_scans_last_before_first_is_zero() {
        assert_eq!(compute_num_scans(10, 3), 0);
    }

    #[test]
    fn compute_num_scans_saturates_instead_of_overflowing() {
        assert_eq!(compute_num_scans(0, u32::MAX), u32::MAX);
        assert_eq!(compute_num_scans(u32::MAX, u32::MAX), 1);
    }

    // Regression test for a second fuzzer-found crash: `data_addr +
    // entry.offset` (both untrusted u64 values) overflowed and panicked for
    // a crafted combination near u64::MAX.
    #[test]
    fn saturating_add_of_addr_and_offset_does_not_panic() {
        let data_addr = u64::MAX - 5;
        let offset = 100u64;
        assert_eq!(data_addr.saturating_add(offset), u64::MAX);
    }

    // Regression test for a third fuzzer-found crash: `read_scan` /
    // `read_scan_flat` / `read_scan_srm_v66` / `read_peaks_only` computed
    // `scan_number - first_scan_number` directly, which underflowed (and
    // panicked under debug assertions) whenever a requested `scan_number`
    // fell below the file's declared `first_scan_number` - reachable simply
    // by asking for an out-of-range scan on a file with a large
    // `first_scan_number`.
    #[test]
    fn scan_number_to_index_normal_case() {
        assert_eq!(scan_number_to_index(5, 1), Some(4));
        assert_eq!(scan_number_to_index(1, 1), Some(0));
    }

    #[test]
    fn scan_number_to_index_below_first_scan_is_none_not_a_panic() {
        assert_eq!(scan_number_to_index(0, 5), None);
        assert_eq!(scan_number_to_index(0, u32::MAX), None);
    }
}
