//! Scan metadata and strict spectrum decoding, independent of XML export.
use std::io::{Read, Seek};

use crate::error::Result;
use crate::scan_event::ScanEvent;
use crate::types::{Activation, MsPower, Polarity};
use crate::RawFileReader;

// -- Structured spectrum record (vendor-neutral, no mzML) --

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
    /// Analyzer used for the precursor scan; needed by mzML CV mapping to
    /// disambiguate CID vs beam-type CID on FTMS instruments.
    pub analyzer: Option<crate::Analyzer>,
}

/// One fully-decoded spectrum with all metadata needed to emit mzML or
/// populate an in-memory record set.
///
/// Returned by [`try_extract_spectrum`] and [`try_iter_spectra`]. Native
/// consumers should use these strict APIs to ingest spectra without an XML
/// serialization round trip. The compatibility helpers [`extract_spectrum`]
/// and [`iter_spectra`] discard decoding errors.
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

/// Extract the metadata of the scan at scan-index `idx` (zero-based offset
/// from the first scan), without reading its peak arrays.
///
/// This is the single place where per-scan fields are derived from the scan
/// index, scan event and trailer. [`extract_spectrum`], the mzML writer, the
/// `openmassspec_core` adapter and the Python bindings all build on it, so a
/// field wired here reaches every output. Returns `None` if `idx` is out of
/// range. Multiplexed reactions remain available in `raw.scan_events`.
/// The singular precursor field is unavailable for these scans.
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

    let precursor = if level == 2 && event.is_none_or(|e| e.reactions.len() <= 1) {
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
                .filter(|&n| n >= first_scan && n < scan_number);
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
/// centroided peak list, matching [`crate::mzml::write_mzml`].
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
    })
}

/// Compatibility extraction that discards decode errors.
/// Use [`try_extract_spectrum`] for conversion or scientific processing.
#[deprecated(
    since = "0.1.0-rc.2",
    note = "Use try_extract_spectrum to retain decoding failures"
)]
pub fn extract_spectrum<R: Read + Seek>(
    raw: &RawFileReader,
    source: &mut R,
    idx: u32,
    include_profile: bool,
) -> Option<SpectrumRecord> {
    try_extract_spectrum(raw, source, idx, include_profile).ok()
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
}

/// Strict arrays with full profile precision, before narrowing to f32 for
/// [`SpectrumRecord`]. Centroid intensities are widened losslessly from f32.
#[derive(Debug, Clone)]
pub struct SpectrumArrays {
    pub metadata: ScanMetadata,
    pub mz: Vec<f64>,
    pub intensity: Vec<f64>,
    pub scan_mode: Option<crate::ScanMode>,
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
    let intensity: Vec<f32> = arrays
        .intensity
        .into_iter()
        .map(|value| value as f32)
        .collect();
    if intensity.iter().any(|value| !value.is_finite()) {
        return Err(
            crate::Error::CorruptData("intensity cannot be represented as finite f32")
                .with_context(decode_context(raw, idx)),
        );
    }
    Ok(DecodedSpectrum {
        metadata: arrays.metadata,
        mz: arrays.mz,
        intensity,
        scan_mode: arrays.scan_mode,
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
    let (mz, intensity, scan_mode) =
        resolve_scan_arrays(raw, source, meta.scan_number, include_profile, event)?;
    if mz.len() != intensity.len()
        || mz.iter().any(|value| !value.is_finite() || *value <= 0.0)
        || intensity.iter().any(|value| !value.is_finite())
    {
        return Err(crate::Error::CorruptData("invalid decoded peak arrays"));
    }
    Ok(SpectrumArrays {
        metadata: meta,
        mz,
        intensity,
        scan_mode,
    })
}

/// Iterate every scan in `raw` as a [`SpectrumRecord`].
///
/// Skipped scans (those for which the peak arrays cannot be decoded) are
/// dropped silently. Use [`try_iter_spectra`] for strict processing. The iterator borrows
/// both `raw` and `source` for its lifetime.
#[deprecated(
    since = "0.1.0-rc.2",
    note = "Use try_iter_spectra to retain decoding failures"
)]
pub fn iter_spectra<'a, R: Read + Seek>(
    raw: &'a RawFileReader,
    source: &'a mut R,
    include_profile: bool,
) -> impl Iterator<Item = SpectrumRecord> + 'a {
    let n = raw.num_scans;
    let mut idx: u32 = 0;
    std::iter::from_fn(move || {
        while idx < n {
            let cur = idx;
            idx += 1;
            if let Ok(rec) = try_extract_spectrum(raw, source, cur, include_profile) {
                return Some(rec);
            }
        }
        None
    })
}
// -- Helpers used by extract_spectrum --

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
            if pairs
                .iter()
                .any(|(mz, intensity)| !mz.is_finite() || *mz <= 0.0 || !intensity.is_finite())
            {
                return Err(crate::Error::CorruptData("invalid decoded profile arrays"));
            }
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
