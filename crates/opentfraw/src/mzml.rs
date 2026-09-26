/// mzML export for Thermo RAW files.
///
/// Writes a valid mzML 1.1.0 document to any `Write` sink. Produces one
/// `<spectrum>` element per scan. Binary arrays (m/z and intensity) are
/// stored as little-endian raw bytes encoded with standard Base64 - no
/// additional compression is applied, keeping this module dependency-free.
///
/// # Usage
/// ```no_run
/// use opentfraw::{RawFileReader, mzml::write_mzml};
/// let raw = RawFileReader::open_path("run.raw").unwrap();
/// let mut out = std::fs::File::create("run.mzML").unwrap();
/// let mut src = std::fs::File::open("run.raw").unwrap();
/// write_mzml(&raw, &mut src, &mut out, "run.raw", false).unwrap();
/// ```
use std::io::{Read, Seek, Write};

use crate::error::Result;
use crate::extra::{scan_extras, ExtraFields};
use crate::types::{Activation, Polarity};
use crate::RawFileReader;

// Preserve the established extraction API paths for downstream callers.
use crate::spectrum::{decode_spectrum, DecodedSpectrum};
pub use crate::spectrum::{
    extract_spectrum, iter_spectra, scan_metadata, try_extract_spectrum, try_iter_spectra,
    PrecursorInfo, ScanMetadata, SpectrumRecord,
};

// -- Adapter / canonical writer wrappers --
//
// The mzML emission machinery itself lives in `openmassspec_core`. Here we
// define a `SpectrumSource` adapter that pulls Thermo scans through
// `extract_spectrum` and converts each opentfraw `SpectrumRecord` into the
// vendor-neutral `openmassspec_core::SpectrumRecord` the canonical writer
// consumes. This keeps opentfraw's public mzML API stable and byte-identical
// to the pre-migration output, while sharing the writer with the other
// vendors.

use openmassspec_core as msc;

const SOFTWARE_NAME: &str = "opentfraw";
const SOFTWARE_VERSION: &str = env!("CARGO_PKG_VERSION");

/// PSI-MS CV term for the source file format (Thermo RAW).
fn source_file_format_cv() -> msc::CvTerm {
    msc::CvTerm::new("MS:1000563", "Thermo RAW format")
}

/// PSI-MS CV term for the native ID format used by Thermo.
fn native_id_format_cv() -> msc::CvTerm {
    msc::CvTerm::new("MS:1000768", "Thermo nativeID format")
}

/// Resolve the instrument CV term for `raw` (lookup mirrors the historical
/// in-crate writer; expands as new Thermo models appear in the corpus).
fn instrument_cv(raw: &RawFileReader) -> msc::CvTerm {
    if let Some(model) = raw.instrument_model {
        let known: &[(&str, &str, &str)] = &[
            ("Orbitrap Astral", "MS:1003378", "Orbitrap Astral"),
            ("Orbitrap Ascend", "MS:1003356", "Orbitrap Ascend"),
            ("Orbitrap Eclipse", "MS:1003029", "Orbitrap Eclipse"),
            (
                "Orbitrap Fusion Lumos",
                "MS:1002732",
                "Orbitrap Fusion Lumos",
            ),
            ("Orbitrap Fusion", "MS:1002416", "Orbitrap Fusion"),
            (
                "Orbitrap Exploris 480",
                "MS:1003028",
                "Orbitrap Exploris 480",
            ),
            (
                "Orbitrap Exploris 240",
                "MS:1003094",
                "Orbitrap Exploris 240",
            ),
            (
                "Orbitrap Exploris 120",
                "MS:1003095",
                "Orbitrap Exploris 120",
            ),
            ("Q Exactive HF-X", "MS:1002877", "Q Exactive HF-X"),
            ("Q Exactive HF", "MS:1002523", "Q Exactive HF"),
            ("Q Exactive Plus", "MS:1002634", "Q Exactive Plus"),
            ("Q Exactive UHMR", "MS:1003245", "Q Exactive UHMR"),
            ("Q Exactive", "MS:1001911", "Q Exactive"),
            ("LTQ Orbitrap Velos Pro", "MS:1003096", "Orbitrap Velos Pro"),
            ("LTQ Orbitrap Velos", "MS:1001742", "LTQ Orbitrap Velos"),
            ("LTQ Orbitrap Elite", "MS:1001910", "Orbitrap Elite"),
            ("LTQ Orbitrap XL", "MS:1000556", "LTQ Orbitrap XL"),
            ("LTQ Orbitrap", "MS:1000449", "LTQ Orbitrap"),
            ("LTQ Velos Pro", "MS:1003495", "Velos Pro"),
            ("LTQ Velos", "MS:1000855", "LTQ Velos"),
            ("LTQ XL", "MS:1000854", "LTQ XL"),
            ("LTQ FT", "MS:1000448", "LTQ FT"),
            ("LTQ", "MS:1000447", "LTQ"),
            ("TSQ Altis", "MS:1002874", "TSQ Altis"),
            ("TSQ Quantiva", "MS:1002418", "TSQ Quantiva"),
            ("TSQ Endura", "MS:1002419", "TSQ Endura"),
            ("TSQ Vantage", "MS:1001510", "TSQ Vantage"),
            ("LCQ Classic", "MS:1000168", "LCQ Classic"),
            ("LCQ Deca", "MS:1000554", "LCQ Deca"),
            ("LCQ Advantage", "MS:1000167", "LCQ Advantage"),
        ];
        for (prefix, acc, name) in known {
            if model.starts_with(prefix) {
                return msc::CvTerm::new(acc, *name);
            }
        }
    }
    msc::CvTerm::new("MS:1000483", "Thermo Fisher Scientific instrument model")
}

/// Acquisition start timestamp (RFC 3339), when the RAW file's info
/// preamble carries a valid date. See
/// [`crate::raw_file_info::RawFileInfoPreamble::acquisition_date_rfc3339`]
/// for the instrument-local-time-with-no-timezone caveat this carries.
fn start_timestamp(raw: &RawFileReader) -> Option<String> {
    raw.raw_file_info.preamble.acquisition_date_rfc3339()
}

fn convert_polarity(p: Option<Polarity>) -> Option<msc::Polarity> {
    p.map(|p| match p {
        Polarity::Negative => msc::Polarity::Negative,
        Polarity::Positive => msc::Polarity::Positive,
    })
}

fn convert_scan_mode(m: Option<crate::ScanMode>) -> Option<msc::ScanMode> {
    m.map(|m| match m {
        crate::ScanMode::Centroid => msc::ScanMode::Centroid,
        crate::ScanMode::Profile => msc::ScanMode::Profile,
    })
}

fn convert_analyzer(a: Option<crate::Analyzer>) -> Option<msc::Analyzer> {
    a.map(|a| match a {
        crate::Analyzer::ITMS => msc::Analyzer::ITMS,
        crate::Analyzer::TQMS => msc::Analyzer::TQMS,
        crate::Analyzer::SQMS => msc::Analyzer::SQMS,
        crate::Analyzer::TOFMS => msc::Analyzer::TOFMS,
        crate::Analyzer::FTMS => msc::Analyzer::FTMS,
        crate::Analyzer::Sector => msc::Analyzer::Sector,
    })
}

fn convert_activation(a: Option<Activation>) -> Option<msc::Activation> {
    a.map(|a| match a {
        Activation::HCD => msc::Activation::HCD,
        Activation::MPID => msc::Activation::MPID,
        Activation::ETD => msc::Activation::ETD,
        Activation::CID => msc::Activation::CID,
        Activation::ECD => msc::Activation::ECD,
        Activation::IRMPD => msc::Activation::IRMPD,
        Activation::PD => msc::Activation::PD,
        Activation::PQD => msc::Activation::PQD,
        Activation::UVPD => msc::Activation::UVPD,
        Activation::SID => msc::Activation::SID,
        Activation::EThcD => msc::Activation::EThcD,
    })
}

/// Thermo native ID string for `scan_number`.
fn native_id_for(scan_number: u32) -> String {
    format!("controllerType=0 controllerNumber=1 scan={scan_number}")
}

fn to_msc_record(
    parts: DecodedSpectrum,
    extra: ::std::collections::BTreeMap<String, String>,
) -> msc::SpectrumRecord {
    let DecodedSpectrum {
        metadata: meta,
        mz,
        intensity,
        scan_mode,
    } = parts;
    let precursor = meta.precursor.map(|p| msc::PrecursorInfo {
        target_mz: p.target_mz,
        selected_mz: p.selected_mz,
        isolation_width: p.isolation_width,
        charge: p.charge,
        intensity: None,
        collision_energy: p.collision_energy,
        ce_is_nce: p.ce_is_nce,
        precursor_native_id: p.master_scan_number.map(native_id_for),
        activation: convert_activation(p.activation),
        analyzer: convert_analyzer(p.analyzer),
        ccs: None,
    });
    msc::SpectrumRecord {
        extra,
        // 0xFFFF is a "no event" sentinel some files store.
        acquisition_event_id: (meta.scan_event != u16::MAX).then_some(u32::from(meta.scan_event)),
        index: meta.index,
        scan_number: meta.scan_number,
        native_id: native_id_for(meta.scan_number),
        ms_level: meta.ms_level,
        polarity: convert_polarity(meta.polarity),
        scan_mode: convert_scan_mode(scan_mode),
        analyzer: convert_analyzer(meta.analyzer),
        filter: meta.filter,
        retention_time_sec: meta.retention_time_min * 60.0,
        total_ion_current: Some(meta.total_ion_current),
        base_peak_mz: Some(meta.base_peak_mz),
        base_peak_intensity: Some(meta.base_peak_intensity),
        low_mz: Some(meta.low_mz),
        high_mz: Some(meta.high_mz),
        ion_injection_time_ms: meta.ion_injection_time_ms,
        inv_mobility: None,
        faims_cv: meta.faims_cv,
        precursor,
        mz,
        intensity,
        inv_mobility_per_peak: None,
    }
}

// -- Chromatogram records (TIC / BPC / SRM) --
//
// These map already-decoded per-scan fields into vendor-neutral
// `openmassspec_core::ChromatogramRecord`s. No new binary parsing happens
// here: the TIC/BPC series come straight from the scan index
// (`total_current`, `base_intensity`), and SRM transitions are formed by
// grouping the same scans by their `scan_event` using the Q1/Q3 maps the
// reader already builds. Following the sibling Waters wiring
// (Sigilweaver/OpenWRaw#9), a trace is only emitted when it can be labeled
// with a real PSI-MS chromatogram-type term and its defining fields resolve;
// anything ambiguous is omitted rather than guessed.

/// Build the total-ion-current chromatogram from `(rt_min, total_current)`
/// pairs, or `None` when the file has no scans.
///
/// `ChromatogramRecord::time_sec` is seconds; the reader reports retention
/// time in minutes, hence the `* 60.0`.
fn tic_record(tic: &[(f64, f64)]) -> Option<msc::ChromatogramRecord> {
    if tic.is_empty() {
        return None;
    }
    Some(msc::ChromatogramRecord {
        index: 0,
        id: "TIC".to_string(),
        chromatogram_type: Some(msc::CvTerm::new(
            "MS:1000235",
            "total ion current chromatogram",
        )),
        precursor_mz: None,
        product_mz: None,
        time_sec: tic.iter().map(|&(rt, _)| (rt * 60.0) as f32).collect(),
        intensity: tic.iter().map(|&(_, i)| i as f32).collect(),
    })
}

/// Build the base-peak chromatogram from `(rt_min, base_intensity, base_mz)`
/// triples, or `None` when the file has no scans. `base_mz` is not carried on
/// a chromatogram record, so only the intensity trace is emitted.
fn bpc_record(bpc: &[(f64, f64, f64)]) -> Option<msc::ChromatogramRecord> {
    if bpc.is_empty() {
        return None;
    }
    Some(msc::ChromatogramRecord {
        index: 0,
        id: "BPC".to_string(),
        chromatogram_type: Some(msc::CvTerm::new("MS:1000628", "basepeak chromatogram")),
        precursor_mz: None,
        product_mz: None,
        time_sec: bpc.iter().map(|&(rt, _, _)| (rt * 60.0) as f32).collect(),
        intensity: bpc.iter().map(|&(_, bpi, _)| bpi as f32).collect(),
    })
}

/// Group SRM/flat-peak scans into one chromatogram per transition.
///
/// `scans` is `(scan_event, rt_min, total_current)` in scan (time) order.
/// Each `scan_event` is one Q1->Q3 transition; its points are the summed
/// product-ion current (`total_current`) already decoded per scan, so no
/// re-summing of peaks is needed. A transition is emitted only when its Q1
/// precursor is known (`q1_by_event`); the product m/z is filled in only when
/// exactly one Q3 window resolves unambiguously, and left unset otherwise
/// rather than guessed. Events are visited in sorted order for deterministic
/// output.
///
/// Two transitions can legitimately share the same Q1/Q3 (e.g. a scheduled
/// method re-monitoring the same pair across separate retention-time
/// windows, each its own `scan_event`), which would otherwise collide on the
/// same `Q1={..} Q3={..}` id - `mzML`'s chromatogram `id` is a plain
/// `xsd:string`, not `xsd:ID`, so nothing upstream would catch a duplicate,
/// but the indexed-mzML writer keys its offset index by `id`
/// (`openmassspec_core::mzml::write_indexed_mzml`), so a collision there
/// would make one of the two transitions unreachable by id-based lookup.
/// Disambiguate any repeat by appending its `scan_event`.
fn srm_chromatograms(
    scans: &[(u16, f64, f64)],
    q1_by_event: &std::collections::HashMap<u16, f64>,
    q3_windows: &std::collections::HashMap<u16, Vec<(f32, f32)>>,
) -> Vec<msc::ChromatogramRecord> {
    let mut points_by_event: std::collections::BTreeMap<u16, Vec<(f32, f32)>> =
        std::collections::BTreeMap::new();
    for &(event, rt_min, current) in scans {
        points_by_event
            .entry(event)
            .or_default()
            .push(((rt_min * 60.0) as f32, current as f32));
    }

    let mut out = Vec::new();
    let mut seen_ids: std::collections::HashSet<String> = std::collections::HashSet::new();
    for (event, points) in points_by_event {
        let Some(&q1) = q1_by_event.get(&event) else {
            continue;
        };
        let product_mz = match q3_windows.get(&event) {
            Some(windows) if windows.len() == 1 => {
                let (lo, hi) = windows[0];
                Some(((lo + hi) / 2.0) as f64)
            }
            _ => None,
        };
        let base_id = match product_mz {
            Some(q3) => format!("SRM Q1={q1:.4} Q3={q3:.4}"),
            None => format!("SRM Q1={q1:.4}"),
        };
        let id = if seen_ids.insert(base_id.clone()) {
            base_id
        } else {
            format!("{base_id} event={event}")
        };
        out.push(msc::ChromatogramRecord {
            index: 0,
            id,
            chromatogram_type: Some(msc::CvTerm::new(
                "MS:1001473",
                "selected reaction monitoring chromatogram",
            )),
            precursor_mz: Some(q1),
            product_mz,
            time_sec: points.iter().map(|&(t, _)| t).collect(),
            intensity: points.iter().map(|&(_, i)| i).collect(),
        });
    }
    out
}

/// Assemble every chromatogram trace `raw` can supply: always a TIC and BPC
/// (when the file has scans), plus one SRM chromatogram per transition for
/// flat-peak SRM files. Records are indexed sequentially in emission order.
fn build_chromatograms(raw: &RawFileReader) -> Vec<msc::ChromatogramRecord> {
    let mut out = Vec::new();
    out.extend(tic_record(&raw.tic_chromatogram()));
    out.extend(bpc_record(&raw.bpc_chromatogram()));

    if raw.flat_peaks && !raw.srm_q1_by_event.is_empty() {
        let scans: Vec<(u16, f64, f64)> = raw
            .scan_index
            .iter()
            .map(|e| (e.scan_event, e.start_time, e.total_current))
            .collect();
        out.extend(srm_chromatograms(
            &scans,
            &raw.srm_q1_by_event,
            &raw.srm_q3_windows,
        ));
    }

    for (i, rec) in out.iter_mut().enumerate() {
        rec.index = i;
    }
    out
}

/// `SpectrumSource` adapter over a Thermo RAW reader.
///
/// Use this when you want to feed Thermo data into any
/// `openmassspec_core`-shaped consumer (the canonical mzML writer, a column
/// store ingester, future Arrow bridge, ...). For the common case of "I just
/// want mzML out", call [`write_mzml`] or [`write_indexed_mzml`] directly.
/// The external `SpectrumSource` trait cannot carry errors. Its iterator stops
/// at the first failure and stores it for [`Self::finish`]. Calling an external
/// canonical writer directly requires checking `finish` afterward.
pub struct OpenTfRawSource<'a, R: Read + Seek> {
    raw: &'a RawFileReader,
    source: &'a mut R,
    raw_filename: &'a str,
    include_profile: bool,
    extra_fields: ExtraFields,
    error: Option<crate::Error>,
}

impl<'a, R: Read + Seek> OpenTfRawSource<'a, R> {
    pub fn new(
        raw: &'a RawFileReader,
        source: &'a mut R,
        raw_filename: &'a str,
        include_profile: bool,
    ) -> Self {
        Self {
            raw,
            source,
            raw_filename,
            include_profile,
            extra_fields: ExtraFields::All,
            error: None,
        }
    }

    /// Return any failure recorded by the infallible vendor-neutral iterator.
    /// Call after fully consuming that iterator, including after a canonical
    /// writer returns success. Prefer this adapter's strict writer methods.
    pub fn finish(&mut self) -> Result<()> {
        self.error.take().map_or(Ok(()), Err)
    }

    /// Write mzML and propagate source failures as well as output failures.
    /// On any error, discard the output because it may contain a partial run.
    pub fn write_mzml<W: Write>(&mut self, out: &mut W) -> Result<()> {
        validate_export(self.raw)?;
        msc::write_mzml(self, out)?;
        self.finish()
    }

    /// Write indexed mzML with strict source error propagation.
    pub fn write_indexed_mzml<W: Write>(&mut self, out: &mut W) -> Result<()> {
        validate_export(self.raw)?;
        msc::write_indexed_mzml(self, out)?;
        self.finish()
    }

    /// Choose which `opentfraw.*` values go into each spectrum's `extra`
    /// map (and so into mzML `<userParam>`s). Defaults to
    /// [`ExtraFields::All`]; see [`crate::extra::extra_field_keys`] for the keys.
    pub fn with_extra_fields(mut self, fields: ExtraFields) -> Self {
        self.extra_fields = fields;
        self
    }
}

/// Distinct analyzers across the file's scan events, in first-seen order.
fn run_analyzers(raw: &RawFileReader) -> Vec<msc::Analyzer> {
    let mut out = Vec::new();
    for event in &raw.scan_events {
        if let Some(a) = convert_analyzer(event.preamble.analyzer()) {
            if !out.contains(&a) {
                out.push(a);
            }
        }
    }
    out
}

impl<'a, R: Read + Seek> msc::SpectrumSource for OpenTfRawSource<'a, R> {
    fn run_metadata(&self) -> msc::RunMetadata {
        msc::RunMetadata {
            extra: ::std::collections::BTreeMap::new(),
            source_file_name: self.raw_filename.to_string(),
            source_file_format: source_file_format_cv(),
            native_id_format: native_id_format_cv(),
            instrument: instrument_cv(self.raw),
            instrument_serial_number: None,
            software_name: SOFTWARE_NAME.into(),
            software_version: SOFTWARE_VERSION.into(),
            acquisition_software_name: None,
            acquisition_software_version: None,
            start_timestamp: start_timestamp(self.raw),
            mobility_array_kind: None,
            analyzers: run_analyzers(self.raw),
        }
    }

    fn iter_spectra<'s>(&'s mut self) -> Box<dyn Iterator<Item = msc::SpectrumRecord> + 's> {
        let n = self.raw.num_scans;
        let raw = self.raw;
        let source = &mut *self.source;
        let include_profile = self.include_profile;
        let extra_fields = &self.extra_fields;
        let error = &mut self.error;
        let mut idx: u32 = 0;
        Box::new(std::iter::from_fn(move || {
            if idx >= n || error.is_some() {
                return None;
            }
            let cur = idx;
            idx += 1;
            match validate_export_scan(raw, cur)
                .and_then(|()| decode_spectrum(raw, source, cur, include_profile))
            {
                Ok(parts) => {
                    let extra = scan_extras(raw, &parts.metadata, extra_fields);
                    Some(to_msc_record(parts, extra))
                }
                Err(err) => {
                    *error = Some(err);
                    None
                }
            }
        }))
    }

    fn spectrum_count_hint(&self) -> Option<usize> {
        Some(self.raw.num_scans as usize)
    }

    fn iter_chromatograms<'s>(
        &'s mut self,
    ) -> Box<dyn Iterator<Item = msc::ChromatogramRecord> + 's> {
        Box::new(build_chromatograms(self.raw).into_iter())
    }
}

fn validate_export_scan(raw: &RawFileReader, idx: u32) -> Result<()> {
    if let Some(event) = raw.scan_events.get(idx as usize) {
        if event.reactions.len() > 1 {
            return Err(crate::Error::UnsupportedOperation(
                "mzML export cannot represent multiple precursor reactions",
            ));
        }
        if event.preamble.bytes.get(6).is_some_and(|&level| level > 2) {
            return Err(crate::Error::UnsupportedOperation(
                "mzML export cannot represent MS3 or higher precursor chains",
            ));
        }
    }
    let scan_number = raw
        .run_header
        .sample_info
        .first_scan_number
        .checked_add(idx)
        .ok_or(crate::Error::InvalidScanIndex(idx))?;
    if raw.scan_parameters(scan_number).is_some_and(|record| {
        record.values.iter().any(|(key, value)| {
            (key.starts_with("SPS Mass ") || key.starts_with("SPS Masses"))
                && match value {
                    crate::generic_data::GenericValue::Float32(v) => *v > 0.0,
                    crate::generic_data::GenericValue::Float64(v) => *v > 0.0,
                    crate::generic_data::GenericValue::String(v) => !v.trim().is_empty(),
                    _ => false,
                }
        })
    }) {
        return Err(crate::Error::UnsupportedOperation(
            "mzML export cannot represent SPS precursors",
        ));
    }
    let meta = scan_metadata(raw, idx).ok_or(crate::Error::InvalidScanIndex(idx))?;
    if meta.ms_level > 1
        && meta
            .precursor
            .as_ref()
            .is_none_or(|p| p.activation.is_none())
    {
        return Err(crate::Error::UnsupportedOperation(
            "mzML writer cannot represent unknown activation without inventing CID",
        ));
    }
    Ok(())
}

fn validate_export(raw: &RawFileReader) -> Result<()> {
    for idx in 0..raw.num_scans {
        validate_export_scan(raw, idx)?;
    }
    Ok(())
}

// -- Public mzML entry points (unchanged signatures) --

/// Write the contents of `raw` as mzML 1.1.0 to `out`.
///
/// * `source` - an open handle to the original `.raw` file (needed to read
///   scan data packets).
/// * `raw_filename` - the file name used for the `<sourceFile>` element.
/// * `include_profile` - when `true`, profile-mode scans export the raw
///   profile signal instead of the centroid peak list.
///
/// All spectra are written. Unreadable peak data returns an error.
/// Discard the output on any error because it may contain a partial run.
pub fn write_mzml<R, W>(
    raw: &RawFileReader,
    source: &mut R,
    out: &mut W,
    raw_filename: &str,
    include_profile: bool,
) -> Result<()>
where
    R: Read + Seek,
    W: Write,
{
    OpenTfRawSource::new(raw, source, raw_filename, include_profile).write_mzml(out)
}

/// Write the contents of `raw` as an indexed mzML 1.1.0 document.
///
/// Indexed mzML adds a `<indexList>` element after all spectra with the byte
/// offset of each `<spectrum>` element, enabling random-access parsing by
/// tools such as pyteomics and pymzml without a full file scan. The
/// `<fileChecksum>` element contains the SHA-1 hash of the file content up
/// to and including `</indexList>`.
pub fn write_indexed_mzml<R, W>(
    raw: &RawFileReader,
    source: &mut R,
    out: &mut W,
    raw_filename: &str,
    include_profile: bool,
) -> Result<()>
where
    R: Read + Seek,
    W: Write,
{
    OpenTfRawSource::new(raw, source, raw_filename, include_profile).write_indexed_mzml(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    #[test]
    fn tic_record_maps_pairs_and_converts_minutes_to_seconds() {
        // (rt_min, total_current)
        let tic = [(0.0, 100.0), (0.5, 250.0), (1.0, 175.0)];
        let rec = tic_record(&tic).expect("non-empty TIC");
        assert_eq!(rec.id, "TIC");
        let cv = rec.chromatogram_type.as_ref().unwrap();
        assert_eq!(cv.accession, "MS:1000235");
        assert_eq!(cv.name, "total ion current chromatogram");
        assert!(rec.precursor_mz.is_none());
        assert!(rec.product_mz.is_none());
        assert_eq!(rec.time_sec, vec![0.0, 30.0, 60.0]); // rt_min * 60
        assert_eq!(rec.intensity, vec![100.0, 250.0, 175.0]);
    }

    #[test]
    fn bpc_record_uses_base_intensity_and_drops_base_mz() {
        // (rt_min, base_intensity, base_mz)
        let bpc = [(0.0, 90.0, 500.0), (0.5, 220.0, 501.0)];
        let rec = bpc_record(&bpc).expect("non-empty BPC");
        assert_eq!(rec.id, "BPC");
        let cv = rec.chromatogram_type.as_ref().unwrap();
        assert_eq!(cv.accession, "MS:1000628");
        assert_eq!(cv.name, "basepeak chromatogram");
        assert_eq!(rec.time_sec, vec![0.0, 30.0]);
        assert_eq!(rec.intensity, vec![90.0, 220.0]); // base_mz not carried
    }

    #[test]
    fn tic_and_bpc_records_are_none_when_no_scans() {
        assert!(tic_record(&[]).is_none());
        assert!(bpc_record(&[]).is_none());
    }

    #[test]
    fn srm_groups_scans_by_event_and_sets_precursor_and_product() {
        // Two transitions interleaved in scan order, one scan_event each.
        let scans = [
            (1u16, 0.0f64, 10.0f64),
            (2u16, 0.0, 20.0),
            (1, 0.5, 12.0),
            (2, 0.5, 22.0),
        ];
        let mut q1 = HashMap::new();
        q1.insert(1u16, 500.0);
        q1.insert(2u16, 600.0);
        let mut q3 = HashMap::new();
        q3.insert(1u16, vec![(100.0f32, 101.0f32)]); // single window -> product m/z
        q3.insert(2u16, vec![(200.0f32, 201.0f32)]);

        let recs = srm_chromatograms(&scans, &q1, &q3);
        assert_eq!(recs.len(), 2);

        // BTreeMap order: event 1 first.
        let r1 = &recs[0];
        assert_eq!(r1.precursor_mz, Some(500.0));
        assert_eq!(r1.product_mz, Some(100.5)); // (100 + 101) / 2
        assert_eq!(r1.id, "SRM Q1=500.0000 Q3=100.5000");
        let cv = r1.chromatogram_type.as_ref().unwrap();
        assert_eq!(cv.accession, "MS:1001473");
        assert_eq!(cv.name, "selected reaction monitoring chromatogram");
        assert_eq!(r1.time_sec, vec![0.0, 30.0]);
        assert_eq!(r1.intensity, vec![10.0, 12.0]);

        let r2 = &recs[1];
        assert_eq!(r2.precursor_mz, Some(600.0));
        assert_eq!(r2.product_mz, Some(200.5));
        assert_eq!(r2.intensity, vec![20.0, 22.0]);
    }

    #[test]
    fn srm_skips_events_without_a_known_q1_and_omits_ambiguous_product() {
        let scans = [
            (1u16, 0.0f64, 10.0f64), // no Q1 -> skipped
            (2u16, 0.0, 20.0),       // Q1 known, two Q3 windows -> product omitted
        ];
        let mut q1 = HashMap::new();
        q1.insert(2u16, 600.0);
        let mut q3 = HashMap::new();
        q3.insert(2u16, vec![(200.0f32, 201.0f32), (300.0f32, 301.0f32)]);

        let recs = srm_chromatograms(&scans, &q1, &q3);
        assert_eq!(recs.len(), 1, "event without a Q1 must be skipped");
        assert_eq!(recs[0].precursor_mz, Some(600.0));
        assert!(
            recs[0].product_mz.is_none(),
            "ambiguous multi-window Q3 must not be guessed"
        );
        assert_eq!(recs[0].id, "SRM Q1=600.0000");
    }

    #[test]
    fn srm_disambiguates_ids_when_two_events_share_the_same_q1_and_q3() {
        // Scheduled method re-monitoring the same transition in two separate
        // time windows: same Q1/Q3, different scan_event.
        let scans = [(1u16, 0.0f64, 10.0f64), (2u16, 5.0, 15.0)];
        let mut q1 = HashMap::new();
        q1.insert(1u16, 500.0);
        q1.insert(2u16, 500.0);
        let mut q3 = HashMap::new();
        q3.insert(1u16, vec![(100.0f32, 101.0f32)]);
        q3.insert(2u16, vec![(100.0f32, 101.0f32)]);

        let recs = srm_chromatograms(&scans, &q1, &q3);
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].id, "SRM Q1=500.0000 Q3=100.5000");
        assert_eq!(
            recs[1].id, "SRM Q1=500.0000 Q3=100.5000 event=2",
            "second event with an identical Q1/Q3 must get a disambiguated id"
        );
        assert_ne!(
            recs[0].id, recs[1].id,
            "ids feed the indexed-mzML offset index and must stay unique"
        );
    }
}
