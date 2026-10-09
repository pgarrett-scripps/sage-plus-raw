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
/// # Example
/// ```no_run
/// use sage_plus_raw::RawFileReader;
/// let raw = RawFileReader::open_path("experiment.raw").unwrap();
/// if let Some(p) = raw.scan_params(1) {
///     println!("Injection time: {:?} ms", p.ion_injection_time_ms());
///     println!("Charge state:   {:?}", p.charge_state());
/// }
/// ```
pub struct ScanParams<'a>(pub &'a GenericRecord);

impl<'a> ScanParams<'a> {
    /// Return the raw `GenericRecord` for direct field access.
    #[inline]
    pub fn record(&self) -> &GenericRecord {
        self.0
    }

    /// Ion injection / fill time in milliseconds.
    ///
    /// Label varies: `"Ion Injection Time (ms):"` (Orbitrap family) vs
    /// `"Ion Inject Time (ms):"` (older LTQ variants).
    pub fn ion_injection_time_ms(&self) -> Option<f64> {
        // Try canonical label first; fall back to legacy label.
        self.0
            .get_f64("Ion Injection Time (ms):")
            .or_else(|| self.0.get_f64("Ion Inject Time (ms):"))
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
        let v = self
            .0
            .get_f64("Monoisotopic M/Z:")
            .or_else(|| self.0.get_f64("MS2 Isolation M/Z:"))
            .or_else(|| self.0.get_f64("Isolation Center M/Z:"))
            .or_else(|| self.0.get_f64("Precursor M/Z:"))?;
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

    /// Number of lock masses found / matched.
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
    /// Label varies: `"MS2 Isolation Width:"` (most common), `"MSn Isolation Width:"`,
    /// or `"Isolation Width (M/Z):"` on some firmware.
    pub fn isolation_width_mz(&self) -> Option<f64> {
        self.0
            .get_f64("MS2 Isolation Width:")
            .or_else(|| self.0.get_f64("MSn Isolation Width:"))
            .or_else(|| self.0.get_f64("Isolation Width (M/Z):"))
            .or_else(|| self.0.get_f64("MS2 Isolation Width (M/Z):"))
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
    /// Tries several label variants present across instrument families.
    /// NCE (normalized collision energy) labels are checked first because
    /// they reflect the user-set method value and are what reference tools
    /// (ThermoRawFileParser, Proteome Discoverer) report.  eV labels are
    /// used as a fallback when no NCE label is present.
    ///
    /// Label priority:
    /// 1. `"HCD Energy:"` / `"HCD Energy V:"` / `"CE:"` - NCE string form
    /// 2. `"Normalized Collision Energy:"` - ion-trap CID NCE
    /// 3. `"HCD Energy (eV):"` - explicit eV label (Q Exactive HF-X, Exploris)
    /// 4. `"HCD Energy eV:"` - eV variant
    /// 5. `"Collision Energy (eV):"` - ITMS CID eV
    pub fn activation_energy(&self) -> Option<f64> {
        // NCE labels: preferred because they match the user-set method value.
        // Skip 0.0 (sentinel for "not set").
        for label in &["HCD Energy:", "HCD Energy V:", "CE:"] {
            if let Some(s) = self.0.get_string(label) {
                if let Ok(v) = s.trim().trim_end_matches('%').parse::<f64>() {
                    if v > 0.0 {
                        return Some(v);
                    }
                }
            }
        }
        if let Some(v) = self
            .0
            .get_f64("Normalized Collision Energy:")
            .filter(|&v| v > 0.0)
        {
            return Some(v);
        }
        // eV labels: used when no NCE label is available.
        if let Some(v) = self.0.get_f64("HCD Energy (eV):").filter(|&v| v > 0.0) {
            return Some(v);
        }
        if let Some(v) = self.0.get_f64("HCD Energy eV:").filter(|&v| v > 0.0) {
            return Some(v);
        }
        self.0
            .get_f64("Collision Energy (eV):")
            .filter(|&v| v > 0.0)
    }

    /// Whether the value returned by [`Self::activation_energy`] is a normalized
    /// collision energy (NCE, dimensionless %) rather than an absolute eV value.
    ///
    /// Returns `true` when `activation_energy` found a value from an NCE label
    /// (`HCD Energy:`, `HCD Energy V:`, `CE:`, or `Normalized Collision Energy:`).
    /// Returns `false` when only eV labels were present or no energy was found.
    pub fn activation_energy_is_nce(&self) -> bool {
        // Returns true if activation_energy() took the NCE path.
        for label in &["HCD Energy:", "HCD Energy V:", "CE:"] {
            if let Some(s) = self.0.get_string(label) {
                if let Ok(v) = s.trim().trim_end_matches('%').parse::<f64>() {
                    if v > 0.0 {
                        return true;
                    }
                }
            }
        }
        self.0
            .get_f64("Normalized Collision Energy:")
            .filter(|&v| v > 0.0)
            .is_some()
    }

    /// Supplemental activation energy for EThcD scans (the HCD component).
    ///
    /// Returns `None` for non-EThcD scans.
    pub fn supplemental_activation_energy(&self) -> Option<f64> {
        if let Some(v) = self.0.get_f64("Supplemental Activation CE:") {
            return Some(v);
        }
        if let Some(s) = self.0.get_string("Supplemental Activation:") {
            return s.trim().trim_end_matches('%').parse::<f64>().ok();
        }
        None
    }

    /// All possible charge states reported by the precursor selection algorithm.
    ///
    /// Returns `None` when the instrument did not report possible charges.
    /// Some firmware stores them as a space-delimited string (e.g. `"2 3"`);
    /// others use a typed integer for the single selected charge.
    pub fn possible_charge_states(&self) -> Option<Vec<u32>> {
        // String variant: "2 3 4"
        if let Some(s) = self.0.get_string("Possible Charge States:") {
            let v: Vec<u32> = s
                .split_whitespace()
                .filter_map(|t| t.parse::<u32>().ok())
                .collect();
            if !v.is_empty() {
                return Some(v);
            }
        }
        // Integer variant (single charge)
        if let Some(c) = self.charge_state() {
            if c > 0 {
                return Some(vec![c as u32]);
            }
        }
        None
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
            GenericValue::String(s) => match s.trim().to_ascii_lowercase().as_str() {
                "on" | "true" => Some(true),
                "off" | "false" => Some(false),
                _ => None,
            },
            _ => None,
        }
    }

    /// S-Lens RF level (V), typically reported on Q Exactive family.
    pub fn s_lens_rf_level(&self) -> Option<f64> {
        self.0.get_f64("S-Lens RF Level:")
    }

    /// AGC fill percentage (0.0-1.0), reported on Q Exactive HF family.
    pub fn agc_fill(&self) -> Option<f64> {
        self.0.get_f64("AGC Fill:")
    }

    /// Orbitrap analyzer temperature (°C), where available.
    pub fn analyzer_temperature(&self) -> Option<f64> {
        self.0.get_f64("Analyzer Temperature:")
    }

    /// PS injection time in milliseconds (pre-scan injection for Q Exactive).
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

    /// Dynamic retention time shift in minutes (Q Exactive HF-X AutoQC).
    pub fn dynamic_rt_shift_min(&self) -> Option<f64> {
        self.0.get_f64("Dynamic RT Shift (min):")
    }

    /// Lock mass correction applied (ppm) - tries several label variants.
    pub fn lock_mass_correction_ppm(&self) -> Option<f64> {
        self.0
            .get_f64("LM Correction (ppm):")
            .or_else(|| self.0.get_f64("LM m/z-Correction (ppm):"))
    }

    /// Number of lock masses found.
    pub fn number_of_lock_masses(&self) -> Option<i32> {
        self.0
            .get_i32("Number of LM Found:")
            .or_else(|| self.0.get_i32("Number of Lock Masses:"))
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

    /// Raw over-fill time T (used for AGC computation).
    pub fn raw_ovft(&self) -> Option<f64> {
        self.0.get_f64("RawOvFtT:")
    }

    /// Error in the isotopic envelope fit (used for charge-state scoring).
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

// -- Status log (instrument log) typed accessor --

/// Typed accessor for a per-scan instrument-status log entry.
///
/// The instrument log records instrument-state values (temperatures, voltages,
/// pressures, etc.) at the time each scan was acquired. The schema varies
/// across instrument models.
pub struct StatusLogEntry<'a>(pub &'a GenericRecord);

impl<'a> StatusLogEntry<'a> {
    /// Return the raw record for direct field access.
    #[inline]
    pub fn record(&self) -> &GenericRecord {
        self.0
    }

    /// Ion injection time in milliseconds (present on Orbitrap family).
    pub fn ion_injection_time_ms(&self) -> Option<f64> {
        self.0
            .get_f64("Ion Injection Time (ms):")
            .or_else(|| self.0.get_f64("Ion Inject Time (ms):"))
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

    /// API (spray) source voltage (V).
    pub fn spray_voltage(&self) -> Option<f64> {
        self.0
            .get_f64("Spray Voltage (V):")
            .or_else(|| self.0.get_f64("Spray Voltage:"))
            .or_else(|| self.0.get_f32("Spray Voltage:").map(f64::from))
    }

    /// Lock mass reference correction (ppm).
    pub fn lock_mass_correction_ppm(&self) -> Option<f64> {
        self.0
            .get_f64("LM Correction (ppm):")
            .or_else(|| self.0.get_f64("LM m/z-Correction (ppm):"))
    }

    /// Capillary temperature (°C).
    pub fn capillary_temperature(&self) -> Option<f64> {
        self.0
            .get_f64("Capillary Temp (°C):")
            .or_else(|| self.0.get_f64("Capillary Temp:"))
            .or_else(|| self.0.get_f32("Capillary Temp:").map(f64::from))
    }

    /// Number of lock masses found.
    pub fn number_of_lock_masses(&self) -> Option<i32> {
        self.0
            .get_i32("Number of LM Found:")
            .or_else(|| self.0.get_i32("Number of Lock Masses:"))
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
