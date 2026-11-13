//! Typed views over trailer and instrument-status records.
use crate::generic_data::{GenericRecord, GenericValue};

// -- High-level typed accessor for scan parameters --

/// Typed accessor for a scan's extra parameters (`ScanParams` stream).
///
/// The underlying [`GenericRecord`] stores named fields whose labels vary
/// slightly across Thermo instrument families. This wrapper normalises the
/// most common labels so callers do not need to hard-code instrument-specific
/// strings.
///
/// Reading these needs a reader opened with the per-scan trailer records
/// ([`crate::OpenOptions::with_records`]); the ten values decoding a spectrum
/// reads are in [`ScanTrailer`] either way.
///
/// # Example
/// ```no_run
/// # fn main() -> Result<(), sage_plus_raw::Error> {
/// use sage_plus_raw::{OpenOptions, RawFileReader};
/// let raw = RawFileReader::open_path_with("experiment.raw", OpenOptions::with_records())?;
/// let p = raw.scan_params(1)?;
/// println!("Injection time: {:?} ms", p.ion_injection_time_ms());
/// println!("Charge state:   {:?}", p.charge_state());
/// # Ok(())
/// # }
/// ```
#[derive(Debug)]
pub struct ScanParams<'a>(pub &'a GenericRecord);

impl<'a> ScanParams<'a> {
    /// Return the raw `GenericRecord` for direct field access.
    #[inline]
    pub fn record(&self) -> &GenericRecord {
        self.0
    }

    /// Ion injection / fill time in milliseconds
    /// (`"Ion Injection Time (ms):"`).
    pub fn ion_injection_time_ms(&self) -> Option<f64> {
        self.0.get_f64("Ion Injection Time (ms):")
    }

    /// Precursor charge state (0 = unknown / MS1 scan).
    pub fn charge_state(&self) -> Option<i32> {
        self.0
            .get_i32("Charge State:")
            // Some LCQ files use UInt8 for charge state.
            .or_else(|| {
                self.0.get("Charge State:").and_then(|v| match v {
                    GenericValue::UInt8(n) => Some(*n as i32),
                    _ => None,
                })
            })
    }

    /// Monoisotopic precursor m/z (0 = not determined).
    ///
    /// Tries multiple label variants for compatibility across instrument families:
    ///
    /// - `"Monoisotopic M/Z:"` - most common (Q Exactive, Orbitrap Fusion)
    /// - `"MS2 Isolation M/Z:"` - some older LTQ firmware
    ///
    /// Returns `None` when the value is absent or zero (not determined).
    pub fn monoisotopic_mz(&self) -> Option<f64> {
        let v = self.0.get_f64("Monoisotopic M/Z:")?;
        if v > 0.0 {
            Some(v)
        } else {
            None
        }
    }

    /// Number of micro-scans averaged into this scan.
    pub fn micro_scan_count(&self) -> Option<i32> {
        self.0.get_i32("Micro Scan Count:")
    }

    /// Scan number of the scan that triggered this dependent scan, from the
    /// `"Master Scan Number:"` trailer. Returns `None` when it is absent.
    ///
    /// `"Master Index:"` is deliberately not used: on Q Exactive and LTQ
    /// Orbitrap files it is a 0/1 flag, not a scan number.
    pub fn master_scan_number(&self) -> Option<i32> {
        self.0.get_i32("Master Scan Number:")
    }

    /// Orbitrap / FT resolving power (e.g. 60000, 120000).
    pub fn ft_resolution(&self) -> Option<i32> {
        self.orbitrap_resolution()
    }

    /// Compatibility alias for [`Self::number_of_lock_masses`], including its
    /// configured-count fallback. Use [`Self::number_of_matched_lock_masses`]
    /// when only actual matches are wanted.
    pub fn number_of_lm_found(&self) -> Option<i32> {
        self.number_of_lock_masses()
    }

    /// Lock-mass m/z correction applied (ppm).
    pub fn lm_correction_ppm(&self) -> Option<f64> {
        self.lock_mass_correction_ppm()
    }

    /// AGC target fill value (ion count).
    pub fn agc_target(&self) -> Option<i32> {
        self.0.get_i32("AGC Target:")
    }

    /// Whether automated gain control (AGC) was active.
    pub fn agc_enabled(&self) -> Option<bool> {
        match self.0.get("AGC:")? {
            GenericValue::Bool(b) => Some(*b),
            GenericValue::String(s) => Some(s.to_ascii_lowercase().contains("on")),
            _ => None,
        }
    }

    /// Elapsed scan time in seconds (Orbitrap instruments only).
    pub fn elapsed_scan_time_s(&self) -> Option<f64> {
        self.0.get_f64("Elapsed Scan Time (sec):")
    }

    /// Maximum allowed ion injection time in milliseconds.
    pub fn max_ion_time_ms(&self) -> Option<f64> {
        self.0.get_f64("Max. Ion Time (ms):")
    }

    /// MSn isolation window width in m/z.
    ///
    /// Read from `"MS2 Isolation Width:"`.
    pub fn isolation_width_mz(&self) -> Option<f64> {
        self.0.get_f64("MS2 Isolation Width:")
    }

    /// MSn isolation window target m/z (the center of the isolation window),
    /// from the `"Target M/Z:"` trailer.
    ///
    /// Most instruments do not write it; callers then fall back to the
    /// event's first reaction `precursor_mz`, which is the target m/z.
    pub fn isolation_target_mz(&self) -> Option<f64> {
        self.0.get_f64("Target M/Z:")
    }

    /// `"MS2 Isolation Offset:"`, a shift in m/z applied to the isolation
    /// window (often 0). It is an offset, not an m/z, so it is never used as
    /// the isolation target.
    pub fn isolation_offset_mz(&self) -> Option<f64> {
        self.0.get_f64("MS2 Isolation Offset:")
    }

    /// Activation energy (eV or %) for the primary activation step.
    ///
    /// NCE (normalized collision energy) labels are checked first, so a scan
    /// that carries both reports its NCE value. The eV label is used only when
    /// no NCE label is present. [`Self::activation_energy_is_nce`] says which
    /// kind was found. Zero is skipped as "not set".
    ///
    /// Label priority:
    /// 1. `"HCD Energy:"` / `"HCD Energy V:"` - NCE, string form
    /// 2. `"HCD Energy eV:"` - eV
    pub fn activation_energy(&self) -> Option<f64> {
        if let Some(v) = self.nce() {
            return Some(v);
        }
        self.0.get_f64("HCD Energy eV:").filter(|&v| v > 0.0)
    }

    /// Whether the value returned by [`Self::activation_energy`] is a normalized
    /// collision energy (NCE, dimensionless %) rather than an absolute eV value.
    ///
    /// Returns `true` when `activation_energy` found a value from an NCE label
    /// (`HCD Energy:` or `HCD Energy V:`). Returns `false` when only the eV
    /// label was present or no energy was found.
    pub fn activation_energy_is_nce(&self) -> bool {
        self.nce().is_some()
    }

    /// First positive NCE value from the string-typed NCE labels. Zero is the
    /// "not set" sentinel and is skipped.
    fn nce(&self) -> Option<f64> {
        ["HCD Energy:", "HCD Energy V:"].iter().find_map(|label| {
            self.0
                .get_string(label)
                .and_then(|s| s.trim().trim_end_matches('%').parse::<f64>().ok())
                .filter(|&v| v > 0.0)
        })
    }

    /// Supplemental activation energy for EThcD scans (the HCD component),
    /// from `"Supplemental Activation CE:"` or `"Supplemental Activation:"`.
    ///
    /// Returns `None` for non-EThcD scans. Neither label occurs in the public
    /// fixtures, and upstream OpenTFRaw b2418f2 removed this accessor for
    /// that reason; it is kept here because the two-clause EThcD filter this
    /// fork renders needs the HCD component of such a scan, and reads the
    /// primary energy only when it is absent.
    pub fn supplemental_activation_energy(&self) -> Option<f64> {
        if let Some(v) = self.0.get_f64("Supplemental Activation CE:") {
            return Some(v);
        }
        if let Some(s) = self.0.get_string("Supplemental Activation:") {
            return s.trim().trim_end_matches('%').parse::<f64>().ok();
        }
        None
    }

    /// Selected precursor charge as a one-element list, or `None` when
    /// [`Self::charge_state`] is absent or not positive.
    ///
    /// `"Possible Charge States:"`, a space-delimited list this read before,
    /// occurs in no corpus file (upstream OpenTFRaw b2418f2).
    pub fn possible_charge_states(&self) -> Option<Vec<u32>> {
        self.charge_state()
            .filter(|&c| c > 0)
            .map(|c| vec![c as u32])
    }

    /// Raw FAIMS compensation voltage in V, including inactive trailer values.
    /// Scientific scan metadata emits it only when FAIMS is explicitly enabled.
    pub fn faims_cv(&self) -> Option<f64> {
        self.0
            .get_f64("FAIMS CV:")
            .or_else(|| self.0.get_f32("FAIMS CV:").map(f64::from))
    }

    /// Whether FAIMS voltage was active for this scan.
    pub fn faims_voltage_on(&self) -> Option<bool> {
        match self.0.get("FAIMS Voltage On:")? {
            GenericValue::Bool(b) => Some(*b),
            GenericValue::String(s) => {
                // Compared in place: this is read once per scan while a file
                // is opened, and lowercasing it allocated a `String` each
                // time.
                let text = s.trim();
                if ["on", "true"].iter().any(|v| text.eq_ignore_ascii_case(v)) {
                    Some(true)
                } else if ["off", "false"]
                    .iter()
                    .any(|v| text.eq_ignore_ascii_case(v))
                {
                    Some(false)
                } else {
                    None
                }
            }
            _ => None,
        }
    }

    /// S-Lens RF level (V), typically reported on Q Exactive family.
    pub fn s_lens_rf_level(&self) -> Option<f64> {
        self.0.get_f64("S-Lens RF Level:")
    }

    /// AGC fill value (`"AGC Fill:"`).
    pub fn agc_fill(&self) -> Option<f64> {
        self.0.get_f64("AGC Fill:")
    }

    /// Orbitrap analyzer temperature (°C), where available.
    pub fn analyzer_temperature(&self) -> Option<f64> {
        self.0.get_f64("Analyzer Temperature:")
    }

    /// `"PS Inj. Time (ms):"` value in milliseconds.
    pub fn ps_injection_time_ms(&self) -> Option<f64> {
        self.0.get_f64("PS Inj. Time (ms):")
    }

    /// Reagent ion injection time in milliseconds (ETD reagent).
    pub fn reagent_ion_injection_time_ms(&self) -> Option<f64> {
        self.0
            .get_f32("Reagent Ion Injection Time (ms):")
            .map(f64::from)
    }

    /// Whether the reagent AGC was active.
    pub fn reagent_ion_agc(&self) -> Option<bool> {
        match self.0.get("Reagent Ion AGC:")? {
            GenericValue::Bool(b) => Some(*b),
            _ => None,
        }
    }

    /// Source CID energy applied in the ion source (eV).
    pub fn source_cid_energy_ev(&self) -> Option<f64> {
        self.0
            .get_f64("Source CID eV:")
            .or_else(|| self.0.get_f32("API Source CID Energy:").map(f64::from))
    }

    /// Dynamic retention time shift in minutes (`"Dynamic RT Shift (min):"`).
    pub fn dynamic_rt_shift_min(&self) -> Option<f64> {
        self.0.get_f64("Dynamic RT Shift (min):")
    }

    /// Lock mass correction applied (ppm) - tries several label variants.
    pub fn lock_mass_correction_ppm(&self) -> Option<f64> {
        self.0
            .get_f64("LM Correction (ppm):")
            .or_else(|| self.0.get_f64("LM m/z-Correction (ppm):"))
    }

    /// Compatibility count: matched lock masses, falling back to the configured
    /// count when the matched-count label is absent or has an unsupported type.
    ///
    /// The two labels describe different quantities. Prefer
    /// [`Self::number_of_matched_lock_masses`] and
    /// [`Self::number_of_configured_lock_masses`] for new code. This accessor
    /// retains its historical fallback behavior.
    pub fn number_of_lock_masses(&self) -> Option<i32> {
        self.0
            .get_i32("Number of LM Found:")
            .or_else(|| self.0.get_i32("Number of Lock Masses:"))
    }

    /// Number of lock-mass peaks found in this record (`Number of LM Found:`).
    ///
    /// Does not fall back to the configured count. Absent, mistyped, or negative
    /// counts return `None`; zero is preserved. Neither zero nor a positive
    /// count establishes whether a correction was applied to this scan.
    pub fn number_of_matched_lock_masses(&self) -> Option<i32> {
        self.0.get_i32("Number of LM Found:").filter(|&n| n >= 0)
    }

    /// Number of configured lock masses (`Number of Lock Masses:`).
    ///
    /// Does not fall back to the matched count. Absent, mistyped, or negative
    /// counts return `None`; zero is preserved. A configured count does not
    /// establish whether any lock masses matched or correction was applied.
    pub fn number_of_configured_lock_masses(&self) -> Option<i32> {
        self.0.get_i32("Number of Lock Masses:").filter(|&n| n >= 0)
    }

    /// Orbitrap resolution setting (not measured, but requested).
    pub fn orbitrap_resolution(&self) -> Option<i32> {
        self.0
            .get_i32("Orbitrap Resolution:")
            .or_else(|| self.0.get_i32("FT Resolution:"))
    }

    /// SPS (Synchronous Precursor Selection) mass for MS3 channel N (0-based index).
    ///
    /// SPS masses are stored as `"SPS Mass 1:"`, `"SPS Mass 2:"`, ... (1-based).
    pub fn sps_mass(&self, channel: usize) -> Option<f32> {
        let label = format!("SPS Mass {}:", channel + 1);
        self.0.get_f32(&label)
    }

    /// Conversion parameter A (Orbitrap m/z conversion polynomial).
    pub fn conversion_parameter_a(&self) -> Option<f64> {
        self.0.get_f64("Conversion Parameter A:")
    }

    /// Conversion parameter B.
    pub fn conversion_parameter_b(&self) -> Option<f64> {
        self.0.get_f64("Conversion Parameter B:")
    }

    /// Conversion parameter C.
    pub fn conversion_parameter_c(&self) -> Option<f64> {
        self.0.get_f64("Conversion Parameter C:")
    }

    /// `"RawOvFtT:"` value, as stored. Its meaning is not decoded.
    pub fn raw_ovft(&self) -> Option<f64> {
        self.0.get_f64("RawOvFtT:")
    }

    /// `"Error in isotopic envelope fit:"` value, as stored.
    pub fn isotopic_fit_error(&self) -> Option<f64> {
        self.0.get_f64("Error in isotopic envelope fit:")
    }

    /// Scan description string (arbitrary text, set by method or real-time software).
    pub fn scan_description(&self) -> Option<&str> {
        self.0.get_string("Scan Description:")
    }

    /// Multi-inject info string (e.g. `"IT=45 "` for ion-trap fill time).
    pub fn multi_inject_info(&self) -> Option<&str> {
        self.0.get_string("Multi Inject Info:")
    }

    /// HCD energy string - raw value as stored (may be `"28.00"`, `"28%"`, or `"N/A"`).
    pub fn hcd_energy(&self) -> Option<&str> {
        self.0
            .get_string("HCD Energy:")
            .or_else(|| self.0.get_string("HCD Energy V:"))
    }
}

// -- Compact per-scan trailer values --

/// The trailer values decoding one scan reads, extracted while the file is
/// opened: what [`crate::scan_metadata`], the rendered scan filter and
/// [`crate::RawFileReader::source_cid_energy_ev`] consult a scan's
/// [`ScanParams`] for.
///
/// A file's trailer stream is one [`GenericRecord`] per scan, each carrying
/// every field of the stream's schema together with that field's label: on a
/// 136k-scan Orbitrap file, 8.6 million `(String, GenericValue)` pairs and
/// about 550 MiB resident while the file is read. This is the part the
/// decoders use, 136 bytes a scan, so a reader keeps only this by default
/// ([`crate::OpenOptions`]) and still decodes byte-identical spectra and
/// filter strings.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
#[non_exhaustive]
pub struct ScanTrailer {
    /// [`ScanParams::isolation_target_mz`].
    pub isolation_target_mz: Option<f64>,
    /// [`ScanParams::monoisotopic_mz`].
    pub monoisotopic_mz: Option<f64>,
    /// [`ScanParams::isolation_width_mz`].
    pub isolation_width_mz: Option<f64>,
    /// [`ScanParams::charge_state`].
    pub charge_state: Option<i32>,
    /// [`ScanParams::activation_energy`].
    pub activation_energy: Option<f64>,
    /// [`ScanParams::activation_energy_is_nce`].
    pub activation_energy_is_nce: bool,
    /// [`ScanParams::supplemental_activation_energy`].
    pub supplemental_activation_energy: Option<f64>,
    /// [`ScanParams::master_scan_number`].
    pub master_scan_number: Option<i32>,
    /// [`ScanParams::ion_injection_time_ms`].
    pub ion_injection_time_ms: Option<f64>,
    /// [`ScanParams::faims_voltage_on`].
    pub faims_voltage_on: Option<bool>,
    /// [`ScanParams::faims_cv`].
    pub faims_cv: Option<f64>,
    /// [`ScanParams::source_cid_energy_ev`].
    pub source_cid_energy_ev: Option<f64>,
    /// Whether the record carried a source-CID energy field at all, however
    /// it was typed: a malformed one is not evidence that the trailer lacks
    /// it, so it still blocks the instrument method's declared value.
    pub source_cid_field_present: bool,
    /// `"Scan Segment:"`, the method event id this scan names.
    pub method_scan_segment: Option<u16>,
    /// `"Scan Event:"`, the method event id this scan names.
    pub method_scan_event: Option<u16>,
}

impl ScanTrailer {
    /// Read one scan's values out of its trailer record.
    pub fn from_params(params: &ScanParams<'_>) -> Self {
        let record = params.record();
        Self {
            isolation_target_mz: params.isolation_target_mz(),
            monoisotopic_mz: params.monoisotopic_mz(),
            isolation_width_mz: params.isolation_width_mz(),
            charge_state: params.charge_state(),
            activation_energy: params.activation_energy(),
            activation_energy_is_nce: params.activation_energy_is_nce(),
            supplemental_activation_energy: params.supplemental_activation_energy(),
            master_scan_number: params.master_scan_number(),
            ion_injection_time_ms: params.ion_injection_time_ms(),
            faims_voltage_on: params.faims_voltage_on(),
            faims_cv: params.faims_cv(),
            source_cid_energy_ev: params.source_cid_energy_ev(),
            source_cid_field_present: record.get("Source CID eV:").is_some()
                || record.get("API Source CID Energy:").is_some(),
            method_scan_segment: crate::instrument_method::trailer_id(record, "Scan Segment:"),
            method_scan_event: crate::instrument_method::trailer_id(record, "Scan Event:"),
        }
    }
}

// -- Status log (instrument log) typed accessor --

/// One instrument status-log record.
#[derive(Debug)]
pub struct StatusLogRecord {
    /// Time the record was written, in minutes (the same scale as scan start
    /// times).
    pub time: f32,
    /// The record's values, labelled by the status-log header.
    pub record: GenericRecord,
}

/// Typed accessor for a per-scan instrument-status log entry.
///
/// The instrument log records instrument-state values (temperatures, voltages,
/// pressures, etc.) over time; a scan sees the last record written at or
/// before its start time. The schema varies across instrument models.
#[derive(Debug)]
pub struct StatusLogEntry<'a>(pub &'a GenericRecord);

impl<'a> StatusLogEntry<'a> {
    /// Return the raw record for direct field access.
    #[inline]
    pub fn record(&self) -> &GenericRecord {
        self.0
    }

    /// Ion injection time in milliseconds (present on Orbitrap family).
    pub fn ion_injection_time_ms(&self) -> Option<f64> {
        self.0.get_f64("Ion Injection Time (ms):")
    }

    /// Orbitrap / FT resolving power setting.
    pub fn ft_resolution(&self) -> Option<i32> {
        self.0
            .get_i32("Orbitrap Resolution:")
            .or_else(|| self.0.get_i32("FT Resolution:"))
    }

    /// FAIMS compensation voltage (V).
    pub fn faims_cv(&self) -> Option<f64> {
        self.0
            .get_f64("FAIMS CV:")
            .or_else(|| self.0.get_f32("FAIMS CV:").map(f64::from))
    }

    /// S-Lens RF level (V).
    pub fn s_lens_rf_level(&self) -> Option<f64> {
        self.0.get_f64("S-Lens RF Level:")
    }

    /// Orbitrap / analyzer temperature (°C).
    pub fn analyzer_temperature(&self) -> Option<f64> {
        self.0
            .get_f64("Analyzer Temperature:")
            .or_else(|| self.0.get_f32("Analyzer Temperature:").map(f64::from))
    }

    /// Lock mass reference correction (ppm).
    pub fn lock_mass_correction_ppm(&self) -> Option<f64> {
        self.0
            .get_f64("LM Correction (ppm):")
            .or_else(|| self.0.get_f64("LM m/z-Correction (ppm):"))
    }

    /// Compatibility count: matched lock masses, falling back to the configured
    /// count when the matched-count label is absent or has an unsupported type.
    ///
    /// The two labels describe different quantities. Prefer
    /// [`Self::number_of_matched_lock_masses`] and
    /// [`Self::number_of_configured_lock_masses`] for new code. This accessor
    /// retains its historical fallback behavior.
    pub fn number_of_lock_masses(&self) -> Option<i32> {
        self.0
            .get_i32("Number of LM Found:")
            .or_else(|| self.0.get_i32("Number of Lock Masses:"))
    }

    /// Number of lock-mass peaks found in this record (`Number of LM Found:`).
    ///
    /// Does not fall back to the configured count. Absent, mistyped, or negative
    /// counts return `None`; zero is preserved. Neither zero nor a positive
    /// count establishes whether a correction was applied to this scan.
    pub fn number_of_matched_lock_masses(&self) -> Option<i32> {
        self.0.get_i32("Number of LM Found:").filter(|&n| n >= 0)
    }

    /// Number of configured lock masses (`Number of Lock Masses:`).
    ///
    /// Does not fall back to the matched count. Absent, mistyped, or negative
    /// counts return `None`; zero is preserved. A configured count does not
    /// establish whether any lock masses matched or correction was applied.
    pub fn number_of_configured_lock_masses(&self) -> Option<i32> {
        self.0.get_i32("Number of Lock Masses:").filter(|&n| n >= 0)
    }

    /// Get any field by name (pass-through to the underlying record).
    pub fn get(&self, label: &str) -> Option<&GenericValue> {
        self.0.get(label)
    }

    /// Get a float64 field by name.
    pub fn get_f64(&self, label: &str) -> Option<f64> {
        self.0.get_f64(label)
    }

    /// Get an int32 field by name.
    pub fn get_i32(&self, label: &str) -> Option<i32> {
        self.0.get_i32(label)
    }

    /// Get a string field by name.
    pub fn get_string(&self, label: &str) -> Option<&str> {
        self.0.get_string(label)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::generic_data::GenericValue;

    /// A record holding exactly the given labels.
    fn record_of(values: &[(&str, GenericValue)]) -> GenericRecord {
        GenericRecord {
            values: values
                .iter()
                .map(|(label, value)| ((*label).into(), value.clone()))
                .collect(),
        }
    }

    /// Every f64 accessor that reads one label, with no fallback label and
    /// no value filtering: (label, accessor).
    #[allow(clippy::type_complexity)]
    const F64_ACCESSORS: &[(&str, fn(&ScanParams) -> Option<f64>)] = &[
        ("Elapsed Scan Time (sec):", |p| p.elapsed_scan_time_s()),
        ("Max. Ion Time (ms):", |p| p.max_ion_time_ms()),
        ("Target M/Z:", |p| p.isolation_target_mz()),
        ("MS2 Isolation Offset:", |p| p.isolation_offset_mz()),
        ("S-Lens RF Level:", |p| p.s_lens_rf_level()),
        ("AGC Fill:", |p| p.agc_fill()),
        ("Analyzer Temperature:", |p| p.analyzer_temperature()),
        ("PS Inj. Time (ms):", |p| p.ps_injection_time_ms()),
        ("Dynamic RT Shift (min):", |p| p.dynamic_rt_shift_min()),
        ("Conversion Parameter A:", |p| p.conversion_parameter_a()),
        ("Conversion Parameter B:", |p| p.conversion_parameter_b()),
        ("Conversion Parameter C:", |p| p.conversion_parameter_c()),
        ("RawOvFtT:", |p| p.raw_ovft()),
        ("Error in isotopic envelope fit:", |p| {
            p.isotopic_fit_error()
        }),
    ];

    /// Every i32 accessor that reads one label, with no fallback.
    #[allow(clippy::type_complexity)]
    const I32_ACCESSORS: &[(&str, fn(&ScanParams) -> Option<i32>)] = &[
        ("Micro Scan Count:", |p| p.micro_scan_count()),
        ("Master Scan Number:", |p| p.master_scan_number()),
        ("AGC Target:", |p| p.agc_target()),
    ];

    #[test]
    fn single_label_accessors_read_their_own_field() {
        for (label, accessor) in F64_ACCESSORS {
            let record = record_of(&[(label, GenericValue::Float64(1.25))]);
            assert_eq!(accessor(&ScanParams(&record)), Some(1.25), "{label}");
            // A field of the wrong type is not read as a number.
            let mistyped = record_of(&[(label, GenericValue::String("1.25".into()))]);
            assert_eq!(accessor(&ScanParams(&mistyped)), None, "{label} mistyped");
            // No other accessor in the table answers for this label.
            for (other_label, other) in F64_ACCESSORS {
                if other_label != label {
                    assert_eq!(
                        other(&ScanParams(&record)),
                        None,
                        "{other_label} vs {label}"
                    );
                }
            }
        }
        for (label, accessor) in I32_ACCESSORS {
            let record = record_of(&[(label, GenericValue::Int32(17))]);
            assert_eq!(accessor(&ScanParams(&record)), Some(17), "{label}");
            let mistyped = record_of(&[(label, GenericValue::Float64(17.0))]);
            assert_eq!(accessor(&ScanParams(&mistyped)), None, "{label} mistyped");
            for (other_label, other) in I32_ACCESSORS {
                if other_label != label {
                    assert_eq!(
                        other(&ScanParams(&record)),
                        None,
                        "{other_label} vs {label}"
                    );
                }
            }
        }
    }

    #[test]
    fn every_accessor_is_none_for_an_empty_record() {
        let record = record_of(&[]);
        let params = ScanParams(&record);
        for (label, accessor) in F64_ACCESSORS {
            assert_eq!(accessor(&params), None, "{label}");
        }
        for (label, accessor) in I32_ACCESSORS {
            assert_eq!(accessor(&params), None, "{label}");
        }
        assert_eq!(params.ion_injection_time_ms(), None);
        assert_eq!(params.charge_state(), None);
        assert_eq!(params.monoisotopic_mz(), None);
        assert_eq!(params.isolation_width_mz(), None);
        assert_eq!(params.ft_resolution(), None);
        assert_eq!(params.orbitrap_resolution(), None);
        assert_eq!(params.lm_correction_ppm(), None);
        assert_eq!(params.lock_mass_correction_ppm(), None);
        assert_eq!(params.number_of_lm_found(), None);
        assert_eq!(params.number_of_lock_masses(), None);
        assert_eq!(params.number_of_matched_lock_masses(), None);
        assert_eq!(params.number_of_configured_lock_masses(), None);
        assert_eq!(params.agc_enabled(), None);
        assert_eq!(params.faims_cv(), None);
        assert_eq!(params.faims_voltage_on(), None);
        assert_eq!(params.reagent_ion_injection_time_ms(), None);
        assert_eq!(params.reagent_ion_agc(), None);
        assert_eq!(params.source_cid_energy_ev(), None);
        assert_eq!(params.activation_energy(), None);
        assert!(!params.activation_energy_is_nce());
        assert_eq!(params.supplemental_activation_energy(), None);
        assert_eq!(params.possible_charge_states(), None);
        assert_eq!(params.sps_mass(0), None);
        assert_eq!(params.scan_description(), None);
        assert_eq!(params.multi_inject_info(), None);
        assert_eq!(params.hcd_energy(), None);
    }

    #[test]
    fn accessors_with_fallback_labels_read_either_spelling() {
        for label in ["Orbitrap Resolution:", "FT Resolution:"] {
            let record = record_of(&[(label, GenericValue::Int32(60_000))]);
            let params = ScanParams(&record);
            assert_eq!(params.orbitrap_resolution(), Some(60_000), "{label}");
            assert_eq!(params.ft_resolution(), Some(60_000), "{label}");
        }
        for label in ["LM Correction (ppm):", "LM m/z-Correction (ppm):"] {
            let record = record_of(&[(label, GenericValue::Float64(-0.25))]);
            let params = ScanParams(&record);
            assert_eq!(params.lock_mass_correction_ppm(), Some(-0.25), "{label}");
            assert_eq!(params.lm_correction_ppm(), Some(-0.25), "{label}");
        }
    }

    /// Each accessor reads the one label the corpus attests, and the variants
    /// that occur in no corpus file are not read at all (upstream OpenTFRaw
    /// b2418f2): a file that does carry one is still read by name through the
    /// trailer record itself.
    #[test]
    fn only_the_corpus_attested_label_of_each_value_is_read() {
        type Read = fn(&ScanParams<'_>) -> Option<f64>;
        let attested: &[(&str, Read)] = &[
            ("Ion Injection Time (ms):", |p| p.ion_injection_time_ms()),
            ("MS2 Isolation Width:", |p| p.isolation_width_mz()),
            ("Monoisotopic M/Z:", |p| p.monoisotopic_mz()),
        ];
        for (label, read) in attested {
            let record = record_of(&[(label, GenericValue::Float64(12.5))]);
            assert_eq!(read(&ScanParams(&record)), Some(12.5), "{label}");
        }
        // Zero still means "not determined" for the precursor m/z.
        let zero = record_of(&[("Monoisotopic M/Z:", GenericValue::Float64(0.0))]);
        assert_eq!(ScanParams(&zero).monoisotopic_mz(), None);

        let dropped: &[(&str, Read)] = &[
            ("Ion Inject Time (ms):", |p| p.ion_injection_time_ms()),
            ("MSn Isolation Width:", |p| p.isolation_width_mz()),
            ("Isolation Width (M/Z):", |p| p.isolation_width_mz()),
            ("MS2 Isolation Width (M/Z):", |p| p.isolation_width_mz()),
            ("MS2 Isolation M/Z:", |p| p.monoisotopic_mz()),
            ("Isolation Center M/Z:", |p| p.monoisotopic_mz()),
            ("Precursor M/Z:", |p| p.monoisotopic_mz()),
            ("Normalized Collision Energy:", |p| p.activation_energy()),
            ("HCD Energy (eV):", |p| p.activation_energy()),
            ("Collision Energy (eV):", |p| p.activation_energy()),
        ];
        for (label, read) in dropped {
            let record = record_of(&[(label, GenericValue::Float64(12.5))]);
            assert_eq!(read(&ScanParams(&record)), None, "{label}");
        }
        // "CE:" was an NCE spelling of its own; it is not read either.
        let ce = record_of(&[("CE:", GenericValue::String("30".into()))]);
        let params = ScanParams(&ce);
        assert_eq!(params.activation_energy(), None);
        assert!(!params.activation_energy_is_nce());
        // Nor is the undecoded charge-state list.
        let listed = record_of(&[(
            "Possible Charge States:",
            GenericValue::String("2 3 4".into()),
        )]);
        assert_eq!(ScanParams(&listed).possible_charge_states(), None);
        // The status log reads the one spelling too.
        let log = record_of(&[("Ion Inject Time (ms):", GenericValue::Float64(5.0))]);
        assert_eq!(StatusLogEntry(&log).ion_injection_time_ms(), None);
        let log = record_of(&[("Ion Injection Time (ms):", GenericValue::Float64(5.0))]);
        assert_eq!(StatusLogEntry(&log).ion_injection_time_ms(), Some(5.0));
    }

    #[test]
    fn charge_state_reads_both_the_int32_and_uint8_encodings() {
        let int32 = record_of(&[("Charge State:", GenericValue::Int32(3))]);
        assert_eq!(ScanParams(&int32).charge_state(), Some(3));
        let uint8 = record_of(&[("Charge State:", GenericValue::UInt8(2))]);
        assert_eq!(ScanParams(&uint8).charge_state(), Some(2));
        // The selected charge is the only charge list there is.
        assert_eq!(ScanParams(&int32).possible_charge_states(), Some(vec![3]));
        assert_eq!(ScanParams(&uint8).possible_charge_states(), Some(vec![2]));
        let zero = record_of(&[("Charge State:", GenericValue::Int32(0))]);
        assert_eq!(ScanParams(&zero).possible_charge_states(), None);
    }

    #[test]
    fn the_lock_mass_counts_keep_their_own_labels() {
        let matched = record_of(&[("Number of LM Found:", GenericValue::Int32(2))]);
        let params = ScanParams(&matched);
        assert_eq!(params.number_of_matched_lock_masses(), Some(2));
        assert_eq!(params.number_of_configured_lock_masses(), None);
        assert_eq!(params.number_of_lock_masses(), Some(2));
        let configured = record_of(&[("Number of Lock Masses:", GenericValue::Int32(3))]);
        let params = ScanParams(&configured);
        assert_eq!(params.number_of_matched_lock_masses(), None);
        assert_eq!(params.number_of_configured_lock_masses(), Some(3));
        // The compatibility accessor falls back to the configured count.
        assert_eq!(params.number_of_lock_masses(), Some(3));
        assert_eq!(params.number_of_lm_found(), Some(3));
        // A negative count is not a count.
        let negative = record_of(&[("Number of LM Found:", GenericValue::Int32(-1))]);
        assert_eq!(ScanParams(&negative).number_of_matched_lock_masses(), None);
    }

    #[test]
    fn the_flag_and_string_valued_accessors_read_their_own_field() {
        let record = record_of(&[
            ("AGC:", GenericValue::Bool(true)),
            ("FAIMS Voltage On:", GenericValue::String("On".into())),
            ("FAIMS CV:", GenericValue::Float32(-45.0)),
            ("Reagent Ion AGC:", GenericValue::Bool(false)),
            (
                "Reagent Ion Injection Time (ms):",
                GenericValue::Float32(7.5),
            ),
            ("Source CID eV:", GenericValue::Float64(12.0)),
            ("SPS Mass 1:", GenericValue::Float32(400.5)),
            ("SPS Mass 2:", GenericValue::Float32(500.5)),
            ("Scan Description:", GenericValue::String("QC".into())),
            ("Multi Inject Info:", GenericValue::String("IT=45 ".into())),
        ]);
        let params = ScanParams(&record);
        assert_eq!(params.agc_enabled(), Some(true));
        assert_eq!(params.faims_voltage_on(), Some(true));
        assert_eq!(params.faims_cv(), Some(-45.0));
        assert_eq!(params.reagent_ion_agc(), Some(false));
        assert_eq!(params.reagent_ion_injection_time_ms(), Some(7.5));
        assert_eq!(params.source_cid_energy_ev(), Some(12.0));
        assert_eq!(params.sps_mass(0), Some(400.5));
        assert_eq!(params.sps_mass(1), Some(500.5));
        assert_eq!(params.sps_mass(2), None);
        assert_eq!(params.scan_description(), Some("QC"));
        assert_eq!(params.multi_inject_info(), Some("IT=45 "));
        // "AGC:" as a string is read for on/off, and an unknown spelling of
        // the FAIMS flag is not guessed at.
        let strings = record_of(&[
            ("AGC:", GenericValue::String("Off".into())),
            ("FAIMS Voltage On:", GenericValue::String("maybe".into())),
        ]);
        let params = ScanParams(&strings);
        assert_eq!(params.agc_enabled(), Some(false));
        assert_eq!(params.faims_voltage_on(), None);
    }

    #[test]
    fn activation_energy_prefers_nce_labels_over_ev_labels() {
        let nce = record_of(&[
            ("HCD Energy:", GenericValue::String("28.00".into())),
            ("HCD Energy eV:", GenericValue::Float64(45.0)),
        ]);
        let params = ScanParams(&nce);
        assert_eq!(params.activation_energy(), Some(28.0));
        assert!(params.activation_energy_is_nce());
        assert_eq!(params.hcd_energy(), Some("28.00"));
        let ev = record_of(&[("HCD Energy eV:", GenericValue::Float64(45.0))]);
        let params = ScanParams(&ev);
        assert_eq!(params.activation_energy(), Some(45.0));
        assert!(!params.activation_energy_is_nce());
        // A percent sign is stripped from either NCE spelling.
        for label in ["HCD Energy:", "HCD Energy V:"] {
            let percent = record_of(&[(label, GenericValue::String("30%".into()))]);
            let params = ScanParams(&percent);
            assert_eq!(params.activation_energy(), Some(30.0), "{label}");
            assert!(params.activation_energy_is_nce(), "{label}");
        }
        // 0.0 is the "not set" sentinel, so the eV label answers instead.
        let zero = record_of(&[
            ("HCD Energy:", GenericValue::String("0.00".into())),
            ("HCD Energy eV:", GenericValue::Float64(35.0)),
        ]);
        let params = ScanParams(&zero);
        assert_eq!(params.activation_energy(), Some(35.0));
        assert!(!params.activation_energy_is_nce());
    }

    #[test]
    fn supplemental_activation_energy_reads_either_encoding() {
        let number = record_of(&[("Supplemental Activation CE:", GenericValue::Float64(20.0))]);
        assert_eq!(
            ScanParams(&number).supplemental_activation_energy(),
            Some(20.0)
        );
        let string = record_of(&[(
            "Supplemental Activation:",
            GenericValue::String("15%".into()),
        )]);
        assert_eq!(
            ScanParams(&string).supplemental_activation_energy(),
            Some(15.0)
        );
    }

    #[test]
    fn master_and_isolation_labels_are_not_misread() {
        let record = GenericRecord {
            values: vec![
                ("Master Index:".into(), GenericValue::Int32(1)),
                ("MS2 Isolation Offset:".into(), GenericValue::Float64(2.5)),
            ],
        };
        let params = ScanParams(&record);
        assert_eq!(params.master_scan_number(), None);
        assert_eq!(params.isolation_target_mz(), None);
        assert_eq!(params.isolation_offset_mz(), Some(2.5));
        let record = GenericRecord {
            values: vec![
                ("Master Scan Number:".into(), GenericValue::Int32(17)),
                ("Target M/Z:".into(), GenericValue::Float64(500.25)),
            ],
        };
        let params = ScanParams(&record);
        assert_eq!(params.master_scan_number(), Some(17));
        assert_eq!(params.isolation_target_mz(), Some(500.25));
    }
}
