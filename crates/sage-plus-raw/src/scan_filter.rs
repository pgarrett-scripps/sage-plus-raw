/// Scan filter string builder - reconstructs a subset of Thermo scan filter
/// syntax for interoperability with downstream tools.
///
/// A Thermo scan filter is a single-line textual summary of a scan's
/// acquisition parameters. It is consumed by virtually every proteomics
/// tool (Proteome Discoverer, MSFragger, MaxQuant, DIA-NN, Skyline,
/// pyteomics, ...) and is the natural key for correlating peptide
/// identifications back to source scans.
///
/// ## Grammar
///
/// Token order follows the "filter line" that the Finnigan Perl module
/// (Gene Selkov, <https://metacpan.org/dist/Finnigan>) renders from
/// `Finnigan::ScanEventPreamble`, `Finnigan::Reaction` and
/// `Finnigan::FractionCollector`.
///
/// ```text
/// <analyzer> <polarity> <scan_mode> <ionization> [<dependent>] <scan_type>
///   ms<n>  [<precursor>@<method><energy> ...]  [<range>]
/// ```
///
/// The standalone builder omits activation. RawFileReader supplies qualified
/// activation when available.
/// These reconstructed filters are not original vendor filter strings.
///
/// An [`Activation::EThcD`] scan renders a single
/// `<precursor>@etd@hcd<supplemental energy>` clause; every other method
/// renders `<precursor>@<method><energy>`.
///
/// ## Examples
///
/// - `FTMS + p NSI Full ms [350.0000-1500.0000]`
/// - `FTMS + c NSI d Full ms2 645.8311 [150.0000-2000.0000]`
/// - `ITMS + c NSI d Full ms2 520.2400 [135.0000-1060.0000]`
/// - `ITMS + c NSI d Full ms2 467.1021@etd130.00 [50.0000-1415.0000]` (ETD)
/// - `FTMS + c NSI d Full ms2 507.8661@etd@hcd25.00 [110.0000-1026.0000]` (EThcD)
/// - `ITMS + c NSI d Full ms3 810.50 265.27 [100.0000-1000.0000]` (MS3)
use std::fmt::Write as _;

use crate::scan_event::ScanEvent;
use crate::scan_index::ScanIndexEntry;
use crate::types::{Activation, Analyzer, MsPower, ScanType};

/// Resolve the activation filter-string token for a given activation code and
/// analyzer. On FTMS instruments, both `CID` (code 4) and `HCD` (code 1) render
/// as "hcd" because they are all beam-type collisions; on ITMS, code 4 renders
/// as "cid".
pub fn activation_str(analyzer: Option<Analyzer>, act: Activation) -> &'static str {
    match act {
        Activation::CID => match analyzer {
            Some(Analyzer::FTMS) => "hcd",
            _ => "cid",
        },
        other => other.as_str(),
    }
}

/// Build the canonical Thermo scan filter string for a single scan event.
///
/// - `event` - the scan event record (provides analyzer, polarity, activation, etc.)
/// - `index_entry` - provides the authoritative m/z scan window
/// - `precursor_mz` - final-stage precursor m/z from scan_params `Monoisotopic M/Z:`
/// - `activation_energy` - primary activation energy (eV or NCE %) from scan_params
/// - `supplemental_energy` - supplemental HCD energy for EThcD scans; `None` for all other types
///
/// For MS2+ scans the function first attempts to build the full precursor chain
/// from `event.reactions` (populated for both pre-v66 and v66 files). If reactions
/// are empty, it falls back to the single `precursor_mz` / `activation_energy` pair.
pub fn build_filter(
    event: &ScanEvent,
    index_entry: &ScanIndexEntry,
    precursor_mz: Option<f64>,
    activation_energy: Option<f64>,
    supplemental_energy: Option<f64>,
) -> String {
    build_filter_with_activation(
        event,
        index_entry,
        precursor_mz,
        activation_energy,
        supplemental_energy,
        None,
        None,
    )
}

/// Build the same filter string as [`build_filter`] into `out`, reusing its
/// allocation instead of returning a new `String`. `out` is cleared first.
///
/// For a loop over many scans this renders every filter into one buffer. The
/// bytes written are exactly those [`build_filter`] returns.
pub fn build_filter_into(
    event: &ScanEvent,
    index_entry: &ScanIndexEntry,
    precursor_mz: Option<f64>,
    activation_energy: Option<f64>,
    supplemental_energy: Option<f64>,
    out: &mut String,
) {
    build_filter_with_activation_into(
        event,
        index_entry,
        precursor_mz,
        activation_energy,
        supplemental_energy,
        None,
        None,
        out,
    )
}

/// Build a filter with a qualified activation and an explicitly sourced
/// source-CID energy. Zero does not establish that source CID is off, so it
/// is omitted rather than rendered as `!sid`.
pub(crate) fn build_filter_with_activation(
    event: &ScanEvent,
    index_entry: &ScanIndexEntry,
    precursor_mz: Option<f64>,
    activation_energy: Option<f64>,
    supplemental_energy: Option<f64>,
    activation: Option<Activation>,
    source_cid_energy_ev: Option<f64>,
) -> String {
    let mut out = String::with_capacity(96);
    build_filter_with_activation_into(
        event,
        index_entry,
        precursor_mz,
        activation_energy,
        supplemental_energy,
        activation,
        source_cid_energy_ev,
        &mut out,
    );
    out
}

/// Render [`build_filter_with_activation`] into `out`, which is cleared
/// first, so a caller rendering many scans keeps one buffer.
///
/// Every field is written straight into `out`: a `write!` into a `String`
/// cannot fail, so the `Result` is discarded rather than unwrapped.
#[allow(clippy::too_many_arguments)]
pub(crate) fn build_filter_with_activation_into(
    event: &ScanEvent,
    index_entry: &ScanIndexEntry,
    precursor_mz: Option<f64>,
    activation_energy: Option<f64>,
    supplemental_energy: Option<f64>,
    activation: Option<Activation>,
    source_cid_energy_ev: Option<f64>,
    out: &mut String,
) {
    out.clear();
    // An unrecognised activation code gets no `@<method>` clause: the raw
    // byte is not a method name (upstream OpenTFRaw eb79dff).
    let activation = activation.filter(|a| !matches!(a, Activation::Unknown(_)));
    let p = &event.preamble;
    let analyzer = p.analyzer();

    // Analyzer (FTMS, ITMS, TQMS, ...)
    if let Some(a) = analyzer {
        out.push_str(a.as_str());
    } else {
        out.push_str("MS");
    }
    out.push(' ');

    // Polarity
    if let Some(pol) = p.polarity() {
        out.push(pol.symbol());
    } else {
        out.push('+');
    }
    out.push(' ');

    // Scan mode (c/p)
    if let Some(m) = p.scan_mode() {
        out.push(m.symbol());
        out.push(' ');
    }

    // Ionization
    if let Some(ion) = p.ionization() {
        out.push_str(ion.as_str());
        out.push(' ');
    }

    // Source CID follows the ionization tokens in the public filter grammar.
    if let Some(energy) = source_cid_energy_ev.filter(|e| e.is_finite() && *e > 0.0) {
        let _ = write!(out, "sid={energy:.2} ");
    }

    // Dependent-scan flag
    if p.is_dependent() {
        out.push_str("d ");
    }

    // Scan type
    match p.scan_type() {
        Some(ScanType::Full) => out.push_str("Full"),
        Some(ScanType::Zoom) => out.push('Z'),
        Some(ScanType::Sim) => out.push_str("SIM"),
        Some(ScanType::Srm) => out.push_str("SRM"),
        Some(ScanType::Crm) => out.push_str("CRM"),
        Some(ScanType::Q1) => out.push_str("Q1MS"),
        Some(ScanType::Q3) => out.push_str("Q3MS"),
        _ => out.push_str("Full"),
    }
    out.push(' ');

    // ms<n> level
    let n = match p.ms_power() {
        Some(MsPower::Ms1) | Some(MsPower::Undefined) | None => 1,
        Some(MsPower::Ms2) => 2,
        Some(MsPower::Ms3) => 3,
        Some(MsPower::Ms4) => 4,
        Some(MsPower::Ms5) => 5,
        Some(MsPower::Ms6) => 6,
        Some(MsPower::Ms7) => 7,
        Some(MsPower::Ms8) => 8,
    };
    if n == 1 {
        out.push_str("ms");
    } else {
        let _ = write!(out, "ms{n}");
    }

    // Preserve reaction targets without inferring energy units or method chains.
    if n >= 2 {
        // Reactions when the event has them, otherwise the single fallback
        // precursor; iterated rather than collected, to keep the render
        // allocation-free.
        let no_reactions = event.reactions.is_empty();
        let masses = event
            .reactions
            .iter()
            .map(|reaction| reaction.precursor_mz)
            .chain(precursor_mz.filter(|_| no_reactions));
        for mz in masses {
            if !mz.is_finite() || mz <= 0.0 {
                continue;
            }
            out.push(' ');
            let _ = write!(out, "{mz:.4}");
            let usable = |value: &f64| value.is_finite() && *value > 0.0;
            match activation {
                // EThcD is written as two activation clauses. The
                // electron-transfer step carries no value of its own and the
                // supplemental HCD energy follows `@hcd`. Filter readers
                // (sage among them) look for `@etd` and `@hcd` separately,
                // so a single `@ethcd` token would match neither.
                Some(Activation::EThcD) => {
                    out.push_str("@etd@hcd");
                    if let Some(energy) = supplemental_energy
                        .filter(usable)
                        .or_else(|| activation_energy.filter(usable))
                    {
                        let _ = write!(out, "{energy:.2}");
                    }
                }
                Some(method) => {
                    out.push('@');
                    out.push_str(method.as_str());
                    if let Some(energy) = activation_energy.filter(usable) {
                        let _ = write!(out, "{energy:.2}");
                    }
                }
                None => {}
            }
        }
    }

    // Scan range [low-high]
    if index_entry.low_mz.is_finite()
        && index_entry.high_mz.is_finite()
        && index_entry.low_mz > 0.0
        && index_entry.high_mz > index_entry.low_mz
    {
        out.push(' ');
        let _ = write!(
            out,
            "[{:.4}-{:.4}]",
            index_entry.low_mz, index_entry.high_mz
        );
    } else if let Some(fc) = event.fraction_collectors.first() {
        if fc.low_mz > 0.0 && fc.high_mz > fc.low_mz {
            out.push(' ');
            let _ = write!(out, "[{:.4}-{:.4}]", fc.low_mz, fc.high_mz);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scan_event::{FractionCollector, Reaction, ScanEvent, ScanEventPreamble};
    use crate::scan_index::ScanIndexEntry;

    fn make_preamble(
        polarity: u8,
        scan_mode: u8,
        ms_power: u8,
        scan_type: u8,
    ) -> ScanEventPreamble {
        let mut bytes = vec![0u8; 136];
        bytes[4] = polarity;
        bytes[5] = scan_mode;
        bytes[6] = ms_power;
        bytes[7] = scan_type;
        bytes[11] = 5; // NSI
        bytes[40] = 4; // FTMS
        ScanEventPreamble { bytes }
    }

    fn make_index(low: f64, high: f64) -> ScanIndexEntry {
        ScanIndexEntry {
            index: 0,
            scan_event: 0,
            scan_segment: 0,
            data_size: 0,
            start_time: 0.0,
            total_current: 0.0,
            base_intensity: 0.0,
            base_mz: 0.0,
            low_mz: low,
            high_mz: high,
            offset: 0,
        }
    }

    #[test]
    fn source_cid_token_placement_and_unknown_values() {
        let ev = ScanEvent {
            preamble: make_preamble(0, 1, 1, 2),
            reactions: vec![],
            fraction_collectors: vec![],
            coefficients: vec![],
        };
        let idx = make_index(100.0, 200.0);
        assert_eq!(
            build_filter_with_activation(&ev, &idx, None, None, None, None, Some(20.0)),
            "FTMS - p NSI sid=20.00 SIM ms [100.0000-200.0000]"
        );
        for energy in [
            None,
            Some(0.0),
            Some(-1.0),
            Some(f64::NAN),
            Some(f64::INFINITY),
        ] {
            assert_eq!(
                build_filter_with_activation(&ev, &idx, None, None, None, None, energy),
                build_filter(&ev, &idx, None, None, None)
            );
        }
    }

    #[test]
    fn ms1_full_profile() {
        let ev = ScanEvent {
            preamble: make_preamble(1, 1, 1, 0),
            reactions: vec![],
            fraction_collectors: vec![FractionCollector {
                low_mz: 350.0,
                high_mz: 1500.0,
            }],
            coefficients: vec![],
        };
        let idx = make_index(350.0, 1500.0);
        let s = build_filter(&ev, &idx, None, None, None);
        assert_eq!(s, "FTMS + p NSI Full ms [350.0000-1500.0000]");
    }

    #[test]
    fn ms2_centroid_hcd_scan_params() {
        let mut pre = make_preamble(1, 0, 2, 0);
        pre.bytes[10] = 1; // dependent
        pre.bytes[24] = 1; // HCD (code 1)
        let ev = ScanEvent {
            preamble: pre,
            reactions: vec![],
            fraction_collectors: vec![FractionCollector {
                low_mz: 150.0,
                high_mz: 2000.0,
            }],
            coefficients: vec![],
        };
        let idx = make_index(150.0, 2000.0);
        let s = build_filter(&ev, &idx, Some(645.8311), Some(28.0), None);
        assert_eq!(s, "FTMS + c NSI d Full ms2 645.8311 [150.0000-2000.0000]");
    }

    #[test]
    fn ms2_reactions_chain() {
        let mut pre = make_preamble(1, 0, 2, 0);
        pre.bytes[10] = 1; // dependent
        pre.bytes[40] = 0; // ITMS
        pre.bytes[24] = 4; // CID
        let rxn = Reaction {
            precursor_mz: 440.254,
            unknown_double: 1.0,
            energy: 35.0,
            unknown_long1: 0,
            unknown_long2: 0,
        };
        let ev = ScanEvent {
            preamble: pre,
            reactions: vec![rxn],
            fraction_collectors: vec![FractionCollector {
                low_mz: 116.0,
                high_mz: 892.0,
            }],
            coefficients: vec![],
        };
        let idx = make_index(116.0, 892.0);
        let s = build_filter(&ev, &idx, Some(440.254), Some(35.0), None);
        // ITMS + CID code 4 → "cid"
        assert_eq!(s, "ITMS + c NSI d Full ms2 440.2540 [116.0000-892.0000]");
    }

    #[test]
    fn ms3_chain() {
        let mut pre = make_preamble(1, 0, 3, 0);
        pre.bytes[10] = 1; // dependent
        pre.bytes[40] = 0; // ITMS
        pre.bytes[24] = 4; // CID
        let ev = ScanEvent {
            preamble: pre,
            reactions: vec![
                Reaction {
                    precursor_mz: 810.50,
                    unknown_double: 1.0,
                    energy: 35.0,
                    unknown_long1: 0,
                    unknown_long2: 0,
                },
                Reaction {
                    precursor_mz: 265.27,
                    unknown_double: 1.0,
                    energy: 35.0,
                    unknown_long1: 0,
                    unknown_long2: 0,
                },
            ],
            fraction_collectors: vec![FractionCollector {
                low_mz: 100.0,
                high_mz: 1000.0,
            }],
            coefficients: vec![],
        };
        let idx = make_index(100.0, 1000.0);
        let s = build_filter(&ev, &idx, Some(265.27), Some(35.0), None);
        assert_eq!(
            s,
            "ITMS + c NSI d Full ms3 810.5000 265.2700 [100.0000-1000.0000]"
        );
    }

    #[test]
    fn ethcd_filter() {
        let mut pre = make_preamble(1, 0, 2, 0);
        pre.bytes[10] = 1; // dependent
        pre.bytes[40] = 4; // FTMS
        pre.bytes[24] = 12; // EThcD
        let ev = ScanEvent {
            preamble: pre,
            reactions: vec![],
            fraction_collectors: vec![FractionCollector {
                low_mz: 150.0,
                high_mz: 2000.0,
            }],
            coefficients: vec![],
        };
        let idx = make_index(150.0, 2000.0);
        let s = build_filter(&ev, &idx, Some(649.12), Some(35.0), Some(28.0));
        assert_eq!(s, "FTMS + c NSI d Full ms2 649.1200 [150.0000-2000.0000]");
    }

    /// v66 MS2 preamble of an electron-transfer event: dependent, a named
    /// code at byte 24, and the electron-transfer flag at byte 120.
    fn etd_preamble(analyzer: u8) -> ScanEventPreamble {
        let mut pre = make_preamble(1, 0, 2, 0);
        pre.bytes[10] = 1; // dependent
        pre.bytes[40] = analyzer;
        pre.bytes[24] = 4; // named code, identical on MS1 and MS2 events
        pre.bytes[120] = 1; // electron transfer
        pre
    }

    fn reaction(precursor_mz: f64, tag: u32) -> Reaction {
        Reaction {
            precursor_mz,
            unknown_double: 0.0,
            energy: 0.0,
            unknown_long1: tag,
            unknown_long2: 0,
        }
    }

    fn event(preamble: ScanEventPreamble, reactions: Vec<Reaction>) -> ScanEvent {
        ScanEvent {
            preamble,
            reactions,
            fraction_collectors: vec![],
            coefficients: vec![],
        }
    }

    #[test]
    fn etd_flag_renders_a_single_etd_clause() {
        let ev = event(etd_preamble(0), vec![reaction(467.1021, 0)]);
        assert_eq!(ev.qualified_activation(66), Some(Activation::ETD));
        let idx = make_index(50.0, 1415.0);
        assert_eq!(
            build_filter_with_activation(
                &ev,
                &idx,
                Some(467.1021),
                Some(130.0),
                None,
                ev.qualified_activation(66),
                None,
            ),
            "ITMS + c NSI d Full ms2 467.1021@etd130.00 [50.0000-1415.0000]"
        );
    }

    #[test]
    fn ethcd_event_renders_two_activation_clauses() {
        let ev = event(
            etd_preamble(4),
            vec![reaction(507.8661, 0), reaction(0.0, 0)],
        );
        assert!(ev.has_supplemental_reaction());
        assert_eq!(ev.qualified_activation(66), Some(Activation::EThcD));
        let idx = make_index(110.0, 1026.0);
        let expected = "FTMS + c NSI d Full ms2 507.8661@etd@hcd25.00 [110.0000-1026.0000]";
        assert_eq!(
            build_filter_with_activation(
                &ev,
                &idx,
                Some(507.8661),
                None,
                Some(25.0),
                ev.qualified_activation(66),
                None,
            ),
            expected
        );
        // Without a supplemental label the primary energy fills the clause:
        // on these events it is the supplemental HCD value.
        assert_eq!(
            build_filter_with_activation(
                &ev,
                &idx,
                Some(507.8661),
                Some(25.0),
                None,
                ev.qualified_activation(66),
                None,
            ),
            expected
        );
    }

    #[test]
    fn etd_flag_needs_a_v66_msn_event_and_a_named_code() {
        // One reaction, so EThcD cannot apply; only the ETD decision varies.
        let reactions = || vec![reaction(467.1021, 0)];
        // Version 63 and 64 preambles are not read for the flag.
        let mut pre = etd_preamble(0);
        pre.bytes.truncate(128);
        assert_eq!(event(pre, reactions()).qualified_activation(66), None);
        // MS1 event.
        let mut pre = etd_preamble(0);
        pre.bytes[6] = 1;
        assert_eq!(event(pre, reactions()).qualified_activation(66), None);
        // Byte 24 zero or out of range, as on misread scan events.
        for code in [0u8, 11, 200] {
            let mut pre = etd_preamble(0);
            pre.bytes[24] = code;
            assert_eq!(event(pre, reactions()).qualified_activation(66), None);
        }
        // Byte 120 = 2 (Q Exactive, Exploris, Astral) is not electron transfer.
        let mut pre = etd_preamble(0);
        pre.bytes[120] = 2;
        assert_eq!(event(pre, reactions()).qualified_activation(66), None);
        // A pre-v66 file is not read for the flag either.
        assert_eq!(
            event(etd_preamble(0), reactions()).qualified_activation(64),
            None
        );
    }

    #[test]
    fn only_the_two_corpus_attested_codes_are_named_methods() {
        // Byte 24 names HCD (1) and CID (4). Every other non-zero code is
        // the raw byte, carries no method name, and renders no clause; the
        // electron-transfer flag is not read for one either.
        assert_eq!(Activation::from_byte(0), None);
        assert_eq!(Activation::from_byte(1), Some(Activation::HCD));
        assert_eq!(Activation::from_byte(4), Some(Activation::CID));
        for code in [2u8, 3, 5, 6, 7, 8, 9, 10, 12, 200, 255] {
            let decoded = Activation::from_byte(code);
            assert_eq!(decoded, Some(Activation::Unknown(code)));
            assert_eq!(decoded.unwrap().as_str(), "unknown");
            assert!(!decoded.unwrap().is_named_code());
            let mut pre = etd_preamble(4);
            pre.bytes[24] = code;
            let ev = event(pre, vec![reaction(467.1021, 0)]);
            assert_eq!(ev.qualified_activation(66), None);
        }
        assert!(Activation::HCD.is_named_code());
        assert!(Activation::CID.is_named_code());
        assert!(!Activation::ETD.is_named_code());
        assert!(!Activation::EThcD.is_named_code());

        // The renderer omits the clause rather than writing the raw byte.
        let idx = make_index(150.0, 2000.0);
        let mut pre = make_preamble(1, 0, 2, 0);
        pre.bytes[10] = 1; // dependent
        pre.bytes[40] = 4; // FTMS
        pre.bytes[24] = 12;
        let ev = ScanEvent {
            preamble: pre,
            reactions: vec![],
            fraction_collectors: vec![],
            coefficients: vec![],
        };
        // This fork reports no activation from the preamble at all, and an
        // `Unknown` reaching the renderer from anywhere else is still no
        // method name.
        assert_eq!(ev.preamble.activation(), None);
        assert_eq!(
            Activation::from_byte(ev.preamble.activation_byte()),
            Some(Activation::Unknown(12))
        );
        for activation in [None, Some(Activation::Unknown(12))] {
            assert_eq!(
                build_filter_with_activation(
                    &ev,
                    &idx,
                    Some(649.12),
                    Some(35.0),
                    Some(28.0),
                    activation,
                    None,
                ),
                "FTMS + c NSI d Full ms2 649.1200 [150.0000-2000.0000]"
            );
        }
    }

    #[test]
    fn the_reaction_tag_outranks_the_electron_transfer_flag() {
        // Tag 11 stays HCD and tag 9 stays ETD whatever byte 120 holds, even
        // in the two-reaction shape that would otherwise read as EThcD.
        for (tag, expected) in [(11u32, Activation::HCD), (9, Activation::ETD)] {
            for flag in [0u8, 1, 2] {
                let mut pre = etd_preamble(4);
                pre.bytes[120] = flag;
                let ev = event(pre, vec![reaction(507.8661, tag), reaction(0.0, tag)]);
                assert_eq!(ev.qualified_activation(66), Some(expected));
            }
        }
    }

    #[test]
    fn cid_and_hcd_filters_are_byte_identical_whatever_byte_120_holds() {
        let idx = make_index(110.0, 1026.0);
        let render = |flag: u8, tag: u32, analyzer: u8| {
            let mut pre = etd_preamble(analyzer);
            pre.bytes[120] = flag;
            let ev = event(pre, vec![reaction(507.8661, tag)]);
            build_filter_with_activation(
                &ev,
                &idx,
                Some(507.8661),
                Some(30.0),
                None,
                ev.qualified_activation(66),
                None,
            )
        };
        // Tag 11 on an Orbitrap: an unchanged `@hcd` clause.
        assert_eq!(
            render(0, 11, 4),
            "FTMS + c NSI d Full ms2 507.8661@hcd30.00 [110.0000-1026.0000]"
        );
        // An ion-trap CID scan the reaction tag leaves unresolved: no clause,
        // and byte 120 = 0 or 2 does not add one.
        for flag in [0u8, 2] {
            assert_eq!(
                render(flag, 7, 0),
                "ITMS + c NSI d Full ms2 507.8661 [110.0000-1026.0000]"
            );
            assert_eq!(render(flag, 11, 4), render(0, 11, 4));
        }
    }

    /// `(event, index entry, precursor, energy, supplemental, activation, sid)`
    /// for one scan of every shape the renderer branches on: an MS1 with only
    /// a collector range, a dependent FTMS MS2 with a qualified HCD clause, an
    /// ITMS MS3 chain, EThcD's two clauses, and a source-CID scan whose
    /// precursor comes from the fallback rather than a reaction.
    #[allow(clippy::type_complexity)]
    fn renderer_cases() -> Vec<(
        ScanEvent,
        ScanIndexEntry,
        Option<f64>,
        Option<f64>,
        Option<f64>,
        Option<Activation>,
        Option<f64>,
    )> {
        let reaction = |precursor_mz: f64| Reaction {
            precursor_mz,
            unknown_double: 1.0,
            energy: 35.0,
            unknown_long1: 0,
            unknown_long2: 0,
        };
        let collector = |low_mz: f64, high_mz: f64| FractionCollector { low_mz, high_mz };

        let ms1 = ScanEvent {
            preamble: make_preamble(1, 1, 1, 0),
            reactions: vec![],
            fraction_collectors: vec![collector(350.0, 1500.0)],
            coefficients: vec![],
        };

        let mut pre = make_preamble(1, 0, 2, 0);
        pre.bytes[10] = 1; // dependent
        pre.bytes[24] = 1; // HCD
        let ms2 = ScanEvent {
            preamble: pre,
            reactions: vec![reaction(645.8311)],
            fraction_collectors: vec![],
            coefficients: vec![],
        };

        let mut pre = make_preamble(1, 0, 3, 0);
        pre.bytes[10] = 1;
        pre.bytes[24] = 4; // CID
        pre.bytes[40] = 0; // ITMS
        let ms3 = ScanEvent {
            preamble: pre,
            reactions: vec![reaction(810.50), reaction(265.27)],
            fraction_collectors: vec![collector(100.0, 1000.0)],
            coefficients: vec![],
        };

        let mut pre = make_preamble(1, 0, 2, 0);
        pre.bytes[10] = 1;
        pre.bytes[24] = 4;
        pre.bytes[120] = 1; // electron transfer
        let ethcd = ScanEvent {
            preamble: pre,
            reactions: vec![reaction(507.8661), reaction(0.0)],
            fraction_collectors: vec![],
            coefficients: vec![],
        };

        let mut pre = make_preamble(0, 0, 2, 2); // negative, SIM
        pre.bytes[11] = 4; // APCI
        let fallback = ScanEvent {
            preamble: pre,
            reactions: vec![],
            fraction_collectors: vec![],
            coefficients: vec![],
        };

        vec![
            (ms1, make_index(350.0, 1500.0), None, None, None, None, None),
            (
                ms2,
                make_index(150.0, 2000.0),
                Some(645.8311),
                Some(28.0),
                None,
                Some(Activation::HCD),
                None,
            ),
            (
                ms3,
                make_index(100.0, 1000.0),
                Some(265.27),
                Some(35.0),
                None,
                Some(Activation::CID),
                None,
            ),
            (
                ethcd,
                make_index(110.0, 1026.0),
                Some(507.8661),
                Some(30.0),
                Some(25.0),
                Some(Activation::EThcD),
                None,
            ),
            (
                fallback,
                make_index(f64::NAN, f64::NAN),
                Some(500.0),
                Some(12.0),
                None,
                None,
                Some(20.0),
            ),
        ]
    }

    #[test]
    fn rendering_into_a_buffer_matches_the_owned_string() {
        for (event, index_entry, precursor, energy, supplemental, activation, sid) in
            renderer_cases()
        {
            let mut out = String::new();
            build_filter_with_activation_into(
                &event,
                &index_entry,
                precursor,
                energy,
                supplemental,
                activation,
                sid,
                &mut out,
            );
            assert_eq!(
                out,
                build_filter_with_activation(
                    &event,
                    &index_entry,
                    precursor,
                    energy,
                    supplemental,
                    activation,
                    sid,
                )
            );

            let mut public = String::new();
            build_filter_into(
                &event,
                &index_entry,
                precursor,
                energy,
                supplemental,
                &mut public,
            );
            assert_eq!(
                public,
                build_filter(&event, &index_entry, precursor, energy, supplemental)
            );
        }
    }

    #[test]
    fn one_buffer_renders_every_scan_without_reallocating() {
        let cases = renderer_cases();
        // Longest first, so every later render fits the capacity it leaves.
        let mut expected: Vec<String> = cases
            .iter()
            .map(
                |(event, index_entry, precursor, energy, supplemental, activation, sid)| {
                    build_filter_with_activation(
                        event,
                        index_entry,
                        *precursor,
                        *energy,
                        *supplemental,
                        *activation,
                        *sid,
                    )
                },
            )
            .collect();
        let longest = expected
            .iter()
            .enumerate()
            .max_by_key(|(_, filter)| filter.len())
            .map(|(i, _)| i)
            .unwrap();
        let mut out = String::with_capacity(expected[longest].len());
        let capacity = out.capacity();
        for (i, (event, index_entry, precursor, energy, supplemental, activation, sid)) in
            cases.iter().enumerate()
        {
            build_filter_with_activation_into(
                event,
                index_entry,
                *precursor,
                *energy,
                *supplemental,
                *activation,
                *sid,
                &mut out,
            );
            // No leftover bytes from the previous, longer scan.
            assert_eq!(out, std::mem::take(&mut expected[i]));
            assert_eq!(out.capacity(), capacity);
        }
    }
}
