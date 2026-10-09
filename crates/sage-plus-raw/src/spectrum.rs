//! Scan metadata and strict spectrum decoding.
use std::io::{Read, Seek};

use crate::error::Result;
use crate::scan_event::ScanEvent;
use crate::types::{Activation, MsPower, Polarity};
use crate::RawFileReader;

/// Precursor metadata for an MS2+ spectrum.
///
/// Built from the per-scan parameter table and/or the scan-event reaction
/// list. `target_mz` is the isolation-window center; `selected_mz` is the
/// monoisotopic-resolved precursor (when available). `collision_energy` is
/// either an absolute eV value or, when `ce_is_nce == true`, a normalized
/// collision energy.
#[derive(Debug, Clone, Default)]
pub struct PrecursorInfo {
    pub target_mz: Option<f64>,
    pub selected_mz: Option<f64>,
    pub isolation_width: Option<f64>,
    pub charge: Option<i32>,
    pub collision_energy: Option<f64>,
    pub ce_is_nce: bool,
    pub master_scan_number: Option<u32>,
    pub activation: Option<Activation>,
    /// Analyzer of this scan, to tell CID from beam-type CID on FTMS
    /// instruments.
    pub analyzer: Option<crate::Analyzer>,
}

/// One fully decoded spectrum with its metadata.
///
/// Returned by [`try_extract_spectrum`] and [`try_iter_spectra`].
#[derive(Debug, Clone)]
pub struct SpectrumRecord {
    pub index: usize,
    pub scan_number: u32,
    pub ms_level: u32,
    pub is_ms1: bool,
    /// Legacy heuristic: true for a nondependent MS2+ scan.
    /// This does not establish DIA and must not classify acquisition mode.
    pub is_dia: bool,
    /// Whether broadband isolation is enabled for this scan.
    pub is_wideband: bool,
    pub polarity: Option<Polarity>,
    /// Effective scan mode after `include_profile` resolution.
    pub scan_mode: Option<crate::ScanMode>,
    pub filter: Option<String>,
    /// Retention time in minutes.
    pub retention_time_min: f64,
    pub total_ion_current: f64,
    pub base_peak_mz: f64,
    pub base_peak_intensity: f64,
    pub low_mz: f64,
    pub high_mz: f64,
    pub ion_injection_time_ms: Option<f64>,
    pub faims_cv: Option<f64>,
    pub precursor: Option<PrecursorInfo>,
    pub mz: Vec<f64>,
    pub intensity: Vec<f32>,
    /// Peaks dropped because their m/z was not finite and positive or their
    /// intensity was not finite or overflowed f32.
    pub dropped_peaks: usize,
}

/// Per-scan metadata: every field of a [`SpectrumRecord`] except the peak
/// arrays, plus the scan-index fields that locate the scan in the method.
///
/// Returned by [`scan_metadata`]. Reading it touches no scan data packet, so
/// it is cheap enough to build for every scan in a file.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct ScanMetadata {
    pub index: usize,
    pub scan_number: u32,
    /// Scan event index, as stored in the scan index. Some files store
    /// `0xFFFF` here.
    pub scan_event: u16,
    /// Scan segment index, as stored in the scan index. Some files store
    /// `0xFFFF` here.
    pub scan_segment: u16,
    /// Size of the scan data packet in bytes, as stored in the scan index.
    pub data_size: u32,
    pub ms_level: u32,
    pub is_ms1: bool,
    /// Legacy heuristic: true for a nondependent MS2+ scan.
    /// This does not establish DIA and must not classify acquisition mode.
    pub is_dia: bool,
    /// Whether broadband isolation is enabled for this scan.
    pub is_wideband: bool,
    pub polarity: Option<Polarity>,
    /// Scan mode as recorded in the scan event.
    pub scan_mode: Option<crate::ScanMode>,
    /// Mass analyzer of this scan.
    pub analyzer: Option<crate::Analyzer>,
    pub filter: Option<String>,
    /// Retention time in minutes.
    pub retention_time_min: f64,
    pub total_ion_current: f64,
    pub base_peak_mz: f64,
    pub base_peak_intensity: f64,
    pub low_mz: f64,
    pub high_mz: f64,
    pub ion_injection_time_ms: Option<f64>,
    pub faims_cv: Option<f64>,
    pub precursor: Option<PrecursorInfo>,
}

fn precursor_masses(
    isolation: Option<f64>,
    monoisotopic: Option<f64>,
    reaction: Option<f64>,
) -> (Option<f64>, Option<f64>) {
    let valid = |x: &f64| x.is_finite() && *x > 0.0;
    let target = isolation.filter(valid).or_else(|| reaction.filter(valid));
    let selected = monoisotopic.filter(valid).or(target);
    (target, selected)
}

/// How far back [`derived_master`] looks for the triggering scan.
const MASTER_SEARCH_SCANS: u32 = 5000;

fn event_level(event: &ScanEvent) -> u32 {
    event.preamble.ms_power().map(ms_level).unwrap_or(1)
}

/// The scan that triggered the level-`level` scan at index `idx` when the
/// trailer names none: the nearest earlier scan one level lower. An MSn
/// scan (n > 2) prefers the nearest MS(n-1) scan of the same cycle whose
/// first reaction isolated the same precursor.
fn derived_master(events: &[ScanEvent], first_scan: u32, idx: u32, level: u32) -> Option<u32> {
    let parent_level = level.checked_sub(1).filter(|&l| l > 0)?;
    let target = (parent_level > 1)
        .then(|| events.get(idx as usize)?.reactions.first())
        .flatten()
        .map(|r| r.precursor_mz);
    let mut nearest = None;
    for j in (idx.saturating_sub(MASTER_SEARCH_SCANS)..idx).rev() {
        let event = events.get(j as usize)?;
        let found = event_level(event);
        if found == parent_level {
            let matches = target.is_none_or(|mz| {
                event
                    .reactions
                    .first()
                    .is_some_and(|r| (r.precursor_mz - mz).abs() <= 0.01)
            });
            if matches {
                return first_scan.checked_add(j);
            }
            nearest = nearest.or(first_scan.checked_add(j));
        } else if found < parent_level {
            break;
        }
    }
    nearest
}

/// Extract the metadata of the scan at scan-index `idx` (zero-based offset
/// from the first scan), without reading its peak arrays.
///
/// This is the single place where per-scan fields are derived from the scan
/// index, scan event and trailer; [`try_extract_spectrum`] builds on it.
/// Returns `None` if `idx` is out of
/// range. Every scan above MS1 has a precursor. For multiplexed (MSX) and
/// multi-stage (EThcD, SPS-MS3) events it describes the first reaction, and
/// every reaction remains available in `raw.scan_events`.
pub fn scan_metadata(raw: &RawFileReader, idx: u32) -> Option<ScanMetadata> {
    if idx >= raw.num_scans {
        return None;
    }
    let first_scan = raw.run_header.sample_info.first_scan_number;
    let scan_number = first_scan.checked_add(idx)?;
    let entry = raw.scan_index.get(idx as usize)?;
    let event = raw.scan_events.get(idx as usize);
    let params = raw.scan_params(scan_number);

    let is_srm = raw.flat_peaks;
    let level = if is_srm {
        2
    } else {
        event
            .and_then(|e| e.preamble.ms_power())
            .map(ms_level)
            .unwrap_or(1)
    };
    let polarity = if is_srm {
        None
    } else {
        event.and_then(|e| e.preamble.polarity())
    };
    let scan_mode = if is_srm {
        Some(crate::ScanMode::Centroid)
    } else {
        event.and_then(|e| e.preamble.scan_mode())
    };
    let filter = raw.scan_filter(scan_number);
    let is_ms1 = !is_srm && level == 1;
    let is_dia = event.is_some_and(|e| e.preamble.is_dia());
    let is_wideband = event.is_some_and(|e| e.preamble.is_wideband());
    let srm_q1 = if is_srm {
        raw.srm_q1_by_event.get(&entry.scan_event).copied()
    } else {
        None
    };
    let srm_ce = if is_srm {
        raw.srm_ce_by_event.get(&entry.scan_event).copied()
    } else {
        None
    };

    let precursor = if level > 1 {
        let info = if let Some(q1) = srm_q1 {
            PrecursorInfo {
                target_mz: Some(q1),
                selected_mz: Some(q1),
                isolation_width: None,
                charge: None,
                collision_energy: srm_ce,
                ce_is_nce: false,
                master_scan_number: None,
                activation: event.and_then(|e| e.qualified_activation(raw.version)),
                analyzer: event.and_then(|e| e.preamble.analyzer()),
            }
        } else {
            let reaction = event.and_then(|e| e.reactions.first());
            let (tm, sm) = precursor_masses(
                params.as_ref().and_then(|p| p.isolation_target_mz()),
                params.as_ref().and_then(|p| p.monoisotopic_mz()),
                reaction.map(|r| r.precursor_mz),
            );
            let iw = params.as_ref().and_then(|p| p.isolation_width_mz());
            let ch = params
                .as_ref()
                .and_then(|p| p.charge_state())
                .filter(|&z| z > 0);
            let ae_from_params = params
                .as_ref()
                .and_then(|p| p.activation_energy())
                .filter(|&e| e > 0.0);
            let ae_is_nce = ae_from_params.is_some()
                && params
                    .as_ref()
                    .map(|p| p.activation_energy_is_nce())
                    .unwrap_or(false);
            // Reaction energy has no qualified unit, so only labeled values are exported.
            let ae = ae_from_params;
            let master = params
                .as_ref()
                .and_then(|p| p.master_scan_number())
                .filter(|&n| n > 0)
                .map(|n| n as u32)
                .filter(|&n| n >= first_scan && n < scan_number)
                .or_else(|| derived_master(&raw.scan_events, first_scan, idx, level));
            PrecursorInfo {
                target_mz: tm,
                selected_mz: sm,
                isolation_width: iw,
                charge: ch,
                collision_energy: ae,
                ce_is_nce: ae_is_nce,
                master_scan_number: master,
                activation: event.and_then(|e| e.qualified_activation(raw.version)),
                analyzer: event.and_then(|e| e.preamble.analyzer()),
            }
        };
        Some(info)
    } else {
        None
    };

    let ion_injection_time_ms = params.as_ref().and_then(|p| p.ion_injection_time_ms());
    let faims_cv = params.as_ref().and_then(|p| {
        (p.faims_voltage_on() == Some(true))
            .then(|| p.faims_cv())
            .flatten()
            .filter(|value| value.is_finite())
    });

    Some(ScanMetadata {
        index: idx as usize,
        scan_number,
        scan_event: entry.scan_event,
        scan_segment: entry.scan_segment,
        data_size: entry.data_size,
        ms_level: level,
        is_ms1,
        is_dia,
        is_wideband,
        polarity,
        scan_mode,
        analyzer: event.and_then(|e| e.preamble.analyzer()),
        filter,
        retention_time_min: entry.start_time,
        total_ion_current: entry.total_current,
        base_peak_mz: entry.base_mz,
        base_peak_intensity: entry.base_intensity,
        low_mz: entry.low_mz,
        high_mz: entry.high_mz,
        ion_injection_time_ms,
        faims_cv,
        precursor,
    })
}

/// Extract a single spectrum's record from `raw` at scan-index `idx`
/// (zero-based offset from the first scan).
///
/// Returns an error if the scan index is invalid or arrays cannot be read.
/// `include_profile` controls
/// whether profile-mode scans return the raw profile signal or the
/// centroided peak list.
pub fn try_extract_spectrum<R: Read + Seek>(
    raw: &RawFileReader,
    source: &mut R,
    idx: u32,
    include_profile: bool,
) -> Result<SpectrumRecord> {
    let DecodedSpectrum {
        metadata: meta,
        mz,
        intensity,
        scan_mode,
        dropped_peaks,
    } = decode_spectrum(raw, source, idx, include_profile)?;
    Ok(SpectrumRecord {
        index: meta.index,
        scan_number: meta.scan_number,
        ms_level: meta.ms_level,
        is_ms1: meta.is_ms1,
        is_dia: meta.is_dia,
        is_wideband: meta.is_wideband,
        polarity: meta.polarity,
        scan_mode,
        filter: meta.filter,
        retention_time_min: meta.retention_time_min,
        total_ion_current: meta.total_ion_current,
        base_peak_mz: meta.base_peak_mz,
        base_peak_intensity: meta.base_peak_intensity,
        low_mz: meta.low_mz,
        high_mz: meta.high_mz,
        ion_injection_time_ms: meta.ion_injection_time_ms,
        faims_cv: meta.faims_cv,
        precursor: meta.precursor,
        mz,
        intensity,
        dropped_peaks,
    })
}

/// Iterate all scans without dropping failures. Each result corresponds to
/// exactly one scan index. Collect into `Result<Vec<_>>` to stop on error.
pub fn try_iter_spectra<'a, R: Read + Seek>(
    raw: &'a RawFileReader,
    source: &'a mut R,
    include_profile: bool,
) -> impl Iterator<Item = Result<SpectrumRecord>> + 'a {
    (0..raw.num_scans).map(move |idx| try_extract_spectrum(raw, source, idx, include_profile))
}

/// A scan's metadata, its peak arrays and the scan mode those arrays are in.
pub(crate) struct DecodedSpectrum {
    pub(crate) metadata: ScanMetadata,
    pub(crate) mz: Vec<f64>,
    pub(crate) intensity: Vec<f32>,
    pub(crate) scan_mode: Option<crate::ScanMode>,
    pub(crate) dropped_peaks: usize,
}

/// Strict arrays with full profile precision, before narrowing to f32 for
/// [`SpectrumRecord`]. Centroid intensities are widened losslessly from f32.
#[derive(Debug, Clone)]
pub struct SpectrumArrays {
    pub metadata: ScanMetadata,
    pub mz: Vec<f64>,
    pub intensity: Vec<f64>,
    pub scan_mode: Option<crate::ScanMode>,
    /// Peaks dropped because their m/z was not finite and positive or their
    /// intensity was not finite or overflowed f32.
    pub dropped_peaks: usize,
}

/// Decode validated arrays without narrowing profile intensities. Metadata's
/// recorded scan mode remains available alongside the effective `scan_mode`.
pub fn try_extract_spectrum_arrays<R: Read + Seek>(
    raw: &RawFileReader,
    source: &mut R,
    idx: u32,
    include_profile: bool,
) -> Result<SpectrumArrays> {
    let result = decode_arrays(raw, source, idx, include_profile);
    result.map_err(|error| error.with_context(decode_context(raw, idx)))
}

fn decode_context(raw: &RawFileReader, idx: u32) -> crate::ErrorContext {
    let number = raw
        .run_header
        .sample_info
        .first_scan_number
        .checked_add(idx);
    let mut context = number
        .map(|number| raw.scan_context(number, crate::Operation::DecodeSpectrum))
        .unwrap_or_else(|| crate::ErrorContext::new(crate::Operation::DecodeSpectrum));
    context.scan_index = Some(idx);
    context
}

pub(crate) fn decode_spectrum<R: Read + Seek>(
    raw: &RawFileReader,
    source: &mut R,
    idx: u32,
    include_profile: bool,
) -> Result<DecodedSpectrum> {
    let arrays = try_extract_spectrum_arrays(raw, source, idx, include_profile)?;
    // Every intensity fits f32: decode_arrays drops peaks that overflow it.
    let intensity: Vec<f32> = arrays
        .intensity
        .into_iter()
        .map(|value| value as f32)
        .collect();
    Ok(DecodedSpectrum {
        metadata: arrays.metadata,
        mz: arrays.mz,
        intensity,
        scan_mode: arrays.scan_mode,
        dropped_peaks: arrays.dropped_peaks,
    })
}

pub(crate) fn validate_scan_event(raw: &RawFileReader, idx: u32) -> Result<()> {
    if !raw.flat_peaks {
        let event = raw
            .scan_events
            .get(idx as usize)
            .ok_or(crate::Error::CorruptData("missing scan event"))?;
        match event.preamble.ms_power() {
            None => {
                return Err(crate::Error::UnsupportedLayout(
                    "unknown scan event MS level",
                ))
            }
            Some(MsPower::Undefined) => {
                return Err(crate::Error::CorruptData("undefined scan event MS level"))
            }
            _ => {}
        }
    }
    Ok(())
}

fn decode_arrays<R: Read + Seek>(
    raw: &RawFileReader,
    source: &mut R,
    idx: u32,
    include_profile: bool,
) -> Result<SpectrumArrays> {
    let meta = scan_metadata(raw, idx).ok_or(crate::Error::InvalidScanIndex(idx))?;
    validate_scan_event(raw, idx)?;
    let event = raw.scan_events.get(idx as usize);
    let (mut mz, mut intensity, scan_mode) =
        resolve_scan_arrays(raw, source, meta.scan_number, include_profile, event)?;
    if mz.len() != intensity.len() {
        return Err(crate::Error::CorruptData("invalid decoded peak arrays"));
    }
    let dropped_peaks = drop_invalid_peaks(&mut mz, &mut intensity);
    Ok(SpectrumArrays {
        metadata: meta,
        mz,
        intensity,
        scan_mode,
        dropped_peaks,
    })
}

fn valid_peak(mz: f64, intensity: f64) -> bool {
    mz.is_finite() && mz > 0.0 && intensity.is_finite() && (intensity as f32).is_finite()
}

/// Remove peaks whose m/z is not finite and positive or whose intensity is
/// not finite or overflows f32, so one bad peak does not cost the scan.
/// Returns how many were removed.
fn drop_invalid_peaks(mz: &mut Vec<f64>, intensity: &mut Vec<f64>) -> usize {
    let before = mz.len();
    let mut keep = 0;
    for i in 0..before {
        if valid_peak(mz[i], intensity[i]) {
            mz[keep] = mz[i];
            intensity[keep] = intensity[i];
            keep += 1;
        }
    }
    mz.truncate(keep);
    intensity.truncate(keep);
    before - keep
}

fn ms_level(power: MsPower) -> u32 {
    match power {
        MsPower::Undefined => 1,
        MsPower::Ms1 => 1,
        MsPower::Ms2 => 2,
        MsPower::Ms3 => 3,
        MsPower::Ms4 => 4,
        MsPower::Ms5 => 5,
        MsPower::Ms6 => 6,
        MsPower::Ms7 => 7,
        MsPower::Ms8 => 8,
    }
}

fn validate_profile_calibration(coeffs: &[f64], analyzer: Option<crate::Analyzer>) -> Result<()> {
    let valid_shape = matches!(coeffs.len(), 4 | 5 | 7)
        || (coeffs.is_empty() && analyzer == Some(crate::Analyzer::ITMS));
    if coeffs.iter().any(|value| !value.is_finite()) {
        return Err(crate::Error::CorruptData("nonfinite profile calibration"));
    }
    if !valid_shape {
        return Err(crate::Error::UnsupportedOperation(
            "missing or unsupported profile calibration",
        ));
    }
    Ok(())
}

/// Resolve the m/z and intensity arrays for a single scan.
///
/// When `include_profile=true` AND the scan packet contains profile data, the
/// profile signal is decoded and returned as the primary arrays (with
/// `effective_scan_mode = Some(ScanMode::Profile)`). Otherwise the centroid
/// peak list is used.
///
/// Returns an error when the scan cannot be read.
fn resolve_scan_arrays<R: Read + Seek>(
    raw: &RawFileReader,
    source: &mut R,
    scan_number: u32,
    include_profile: bool,
    event: Option<&ScanEvent>,
) -> Result<(Vec<f64>, Vec<f64>, Option<crate::ScanMode>)> {
    if include_profile && !raw.flat_peaks {
        let packet = raw.read_scan(source, scan_number)?;
        if let Some(profile) = packet.profile {
            let coeffs = event.map(|e| e.coefficients.as_slice()).unwrap_or(&[]);
            validate_profile_calibration(coeffs, event.and_then(|e| e.preamble.analyzer()))?;
            let pairs = profile.to_mz_intensity(coeffs);
            let (mz, int) = pairs.into_iter().unzip();
            return Ok((mz, int, Some(crate::ScanMode::Profile)));
        }
        let mz: Vec<f64> = packet.peaks.iter().map(|p| p.mz).collect();
        let int: Vec<f64> = packet
            .peaks
            .iter()
            .map(|p| f64::from(p.abundance))
            .collect();
        return Ok((mz, int, Some(crate::ScanMode::Centroid)));
    }
    let peaks = raw.read_peaks_only(source, scan_number)?;
    let mz: Vec<f64> = peaks.iter().map(|p| p.mz).collect();
    let int: Vec<f64> = peaks.iter().map(|p| f64::from(p.abundance)).collect();
    Ok((mz, int, Some(crate::ScanMode::Centroid)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn invalid_peaks_are_dropped_not_the_scan() {
        let mut mz = vec![100.0, -1.0, 0.0, f64::NAN, 200.0, 300.0, f64::INFINITY];
        let mut intensity = vec![1.0, 2.0, 3.0, 4.0, f64::NAN, 1e300, 5.0];
        assert_eq!(drop_invalid_peaks(&mut mz, &mut intensity), 6);
        assert_eq!(mz, [100.0]);
        assert_eq!(intensity, [1.0]);
        let mut mz = vec![100.0, 200.0];
        let mut intensity = vec![0.0, 1.0];
        assert_eq!(drop_invalid_peaks(&mut mz, &mut intensity), 0);
        assert_eq!(mz.len(), 2);
    }

    #[test]
    fn profile_calibration_rejects_missing_ftms_and_nonfinite_values() {
        assert!(validate_profile_calibration(&[], Some(crate::Analyzer::FTMS)).is_err());
        assert!(validate_profile_calibration(&[], None).is_err());
        assert!(validate_profile_calibration(&[], Some(crate::Analyzer::ITMS)).is_ok());
        assert!(validate_profile_calibration(&[0.0, 1.0], Some(crate::Analyzer::FTMS)).is_err());
        assert!(validate_profile_calibration(
            &[0.0, 0.0, f64::NAN, 1.0],
            Some(crate::Analyzer::FTMS)
        )
        .is_err());
    }

    fn event(level: u8, mz: Option<f64>) -> ScanEvent {
        let mut bytes = vec![0u8; 136];
        bytes[6] = level;
        ScanEvent {
            preamble: crate::scan_event::ScanEventPreamble { bytes },
            reactions: mz
                .into_iter()
                .map(|precursor_mz| crate::scan_event::Reaction {
                    precursor_mz,
                    unknown_double: 0.0,
                    energy: 0.0,
                    unknown_long1: 0,
                    unknown_long2: 0,
                })
                .collect(),
            fraction_collectors: vec![],
            coefficients: vec![],
        }
    }

    #[test]
    fn derived_master_follows_the_cycle() {
        // MS1, MS2 500, MS2 600, MS3 of 500, MS3 of 600, MS1, MS2 700
        let events = [
            event(1, None),
            event(2, Some(500.0)),
            event(2, Some(600.0)),
            event(3, Some(500.0)),
            event(3, Some(600.0)),
            event(1, None),
            event(2, Some(700.0)),
            event(3, Some(900.0)),
        ];
        let master = |idx| derived_master(&events, 10, idx, event_level(&events[idx as usize]));
        assert_eq!(master(0), None);
        assert_eq!(master(1), Some(10));
        assert_eq!(master(2), Some(10));
        assert_eq!(master(3), Some(11));
        assert_eq!(master(4), Some(12));
        assert_eq!(master(6), Some(15));
        // No same-precursor MS2 in the cycle: the nearest one.
        assert_eq!(master(7), Some(16));
    }

    #[test]
    fn isolation_target_is_not_inferred_from_monoisotopic_mass() {
        assert_eq!(
            precursor_masses(None, Some(499.8), Some(500.0)),
            (Some(500.0), Some(499.8))
        );
        assert_eq!(
            precursor_masses(None, Some(499.8), None),
            (None, Some(499.8))
        );
        assert_eq!(
            precursor_masses(Some(501.0), Some(f64::NAN), Some(500.0)),
            (Some(501.0), Some(501.0))
        );
    }
}
