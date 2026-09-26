//! Independent ProteoWizard acquisition regressions, required by CI.
use openmassspec_core::SpectrumSource;
use opentfraw::{mzml::OpenTfRawSource, Activation, RawFileReader, ScanMode};
use std::{fs::File, io::Cursor, path::PathBuf};

fn fixture(name: &str) -> PathBuf {
    let root =
        PathBuf::from(std::env::var("THERMO_PUBLIC_CORPUS").expect("set THERMO_PUBLIC_CORPUS"));
    let root = if root.is_absolute() {
        root
    } else {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../..")
            .join(root)
    };
    root.join(name)
}

#[test]
#[ignore = "requires downloaded public fixtures, run by CI"]
fn qualified_activation_uses_reaction_tags_not_identical_preambles() {
    let etd = RawFileReader::open_path(fixture("BSA-FT-ETD.raw")).unwrap();
    let hcd = RawFileReader::open_path(fixture("BSA-FT-HCD.raw")).unwrap();
    assert_eq!(
        etd.scan_events[0].preamble.bytes,
        hcd.scan_events[0].preamble.bytes
    );
    for (raw, name, activation, cv, token) in [
        (&etd, "BSA-FT-ETD", Activation::ETD, "MS:1000598", "@etd"),
        (&hcd, "BSA-FT-HCD", Activation::HCD, "MS:1000422", "@hcd"),
    ] {
        let reference = std::fs::read_to_string(fixture(&format!("{name}.mzML"))).unwrap();
        assert!(reference.contains(cv));
        assert!(reference.contains(token));
        let meta = opentfraw::scan_metadata(raw, 0).unwrap();
        assert_eq!(
            meta.precursor.as_ref().unwrap().activation,
            Some(activation)
        );
        assert!(meta.precursor.unwrap().master_scan_number.is_none());
        assert!(meta.filter.unwrap().contains(token));
        let mut output = Vec::new();
        opentfraw::write_mzml(
            raw,
            &mut File::open(fixture(&format!("{name}.raw"))).unwrap(),
            &mut output,
            name,
            false,
        )
        .unwrap();
        let xml = String::from_utf8(output).unwrap();
        assert!(xml.contains(cv));
        assert!(!xml.contains("accession=\"MS:1000133\""));
    }
}

#[test]
#[ignore = "requires downloaded public fixtures, run by CI"]
fn strict_paths_propagate_decode_errors_and_label_centroids() {
    let path = fixture("BSA-FT-HCD.raw");
    let raw = RawFileReader::open_path(&path).unwrap();
    assert_eq!(
        opentfraw::scan_metadata(&raw, 0).unwrap().scan_mode,
        Some(ScanMode::Profile)
    );
    let spectrum =
        opentfraw::try_extract_spectrum(&raw, &mut File::open(&path).unwrap(), 0, false).unwrap();
    assert_eq!(spectrum.scan_mode, Some(ScanMode::Centroid));
    assert!(!spectrum.mz.is_empty());
    assert!(
        opentfraw::try_extract_spectrum(&raw, &mut Cursor::new(Vec::<u8>::new()), 0, false)
            .is_err()
    );
    assert!(
        opentfraw::try_extract_spectrum(&raw, &mut File::open(&path).unwrap(), 1, false).is_err()
    );
    let results: Vec<_> =
        opentfraw::try_iter_spectra(&raw, &mut Cursor::new(Vec::<u8>::new()), false).collect();
    assert_eq!(results.len(), raw.num_scans as usize);
    assert!(results.iter().all(Result::is_err));
    let mut output = Vec::new();
    assert!(opentfraw::write_mzml(
        &raw,
        &mut Cursor::new(Vec::<u8>::new()),
        &mut output,
        "broken.raw",
        false
    )
    .is_err());
    output.clear();
    assert!(opentfraw::write_indexed_mzml(
        &raw,
        &mut Cursor::new(Vec::<u8>::new()),
        &mut output,
        "broken.raw",
        false
    )
    .is_err());
    let mut empty = Cursor::new(Vec::<u8>::new());
    let mut adapter = OpenTfRawSource::new(&raw, &mut empty, "broken.raw", false);
    assert_eq!(adapter.iter_spectra().count(), 0);
    assert!(adapter.finish().is_err());
}

#[test]
#[ignore = "requires downloaded public fixtures, run by CI"]
fn unknown_activation_and_complex_precursors_cannot_be_exported_lossily() {
    let path = fixture("BSA-FT-HCD.raw");
    let mut raw = RawFileReader::open_path(&path).unwrap();
    raw.scan_events[0].reactions[0].unknown_long1 = 999;
    assert!(opentfraw::scan_metadata(&raw, 0)
        .unwrap()
        .precursor
        .unwrap()
        .activation
        .is_none());
    assert!(!raw.scan_filter(1).unwrap().contains('@'));
    let mut output = Vec::new();
    assert!(opentfraw::write_mzml(
        &raw,
        &mut File::open(&path).unwrap(),
        &mut output,
        "unknown.raw",
        false
    )
    .is_err());
    assert!(output.is_empty());
    raw.scan_events[0].reactions[0].unknown_long1 = 11;
    raw.scan_events[0].preamble.bytes[6] = 3;
    assert!(opentfraw::scan_metadata(&raw, 0)
        .unwrap()
        .precursor
        .is_none());
    assert!(opentfraw::write_indexed_mzml(
        &raw,
        &mut File::open(&path).unwrap(),
        &mut output,
        "ms3.raw",
        false
    )
    .is_err());
    assert!(output.is_empty());
    raw.scan_events[0].preamble.bytes[6] = 2;
    raw.scan_parameters[0].values.push((
        "SPS Masses:".into(),
        opentfraw::generic_data::GenericValue::String("677.31".into()),
    ));
    assert!(opentfraw::write_mzml(
        &raw,
        &mut File::open(&path).unwrap(),
        &mut output,
        "sps.raw",
        false
    )
    .is_err());
    assert!(output.is_empty());
}

#[test]
#[ignore = "requires downloaded public fixtures, run by CI"]
fn profile_requires_calibration_and_strict_metadata_requires_an_event() {
    let path = fixture("BSA-FT-HCD.raw");
    let mut raw = RawFileReader::open_path(&path).unwrap();
    let profile =
        opentfraw::try_extract_spectrum(&raw, &mut File::open(&path).unwrap(), 0, true).unwrap();
    assert_eq!(profile.scan_mode, Some(ScanMode::Profile));
    assert!(!profile.mz.is_empty());
    raw.scan_events[0].coefficients.clear();
    assert!(
        opentfraw::try_extract_spectrum(&raw, &mut File::open(&path).unwrap(), 0, true).is_err()
    );
    raw.scan_events.clear();
    assert!(
        opentfraw::try_extract_spectrum(&raw, &mut File::open(&path).unwrap(), 0, false).is_err()
    );
}

#[test]
#[ignore = "requires downloaded public fixtures, run by CI"]
fn multiple_reactions_and_offset_isolation_are_preserved() {
    let path = fixture("FT-HCD-MSX.raw");
    let raw = RawFileReader::open_path(&path).unwrap();
    let expected = [522.987607, 698.067188, 502.978512, 768.099021];
    assert_eq!(raw.scan_events[0].reactions.len(), expected.len());
    for (reaction, expected) in raw.scan_events[0].reactions.iter().zip(expected) {
        assert!((reaction.precursor_mz - expected).abs() < 1e-9);
        assert_eq!(reaction.unknown_double, 5.0);
    }
    assert!(opentfraw::scan_metadata(&raw, 0)
        .unwrap()
        .precursor
        .is_none());
    let mut output = Vec::new();
    assert!(opentfraw::write_mzml(
        &raw,
        &mut File::open(path).unwrap(),
        &mut output,
        "msx.raw",
        false
    )
    .is_err());
    assert!(output.is_empty());
    let raw =
        RawFileReader::open_path(fixture("IsolationMzOffset-ReportedMassOffset.raw")).unwrap();
    let precursor = opentfraw::scan_metadata(&raw, 0)
        .unwrap()
        .precursor
        .unwrap();
    assert!((precursor.target_mz.unwrap() - 379.7113342285156).abs() < 1e-9);
    assert!((precursor.selected_mz.unwrap() - 384.7153374006682).abs() < 1e-9);
    for name in ["IT-HCD-SPS.raw", "source_cid_test_3scans.raw"] {
        assert!(RawFileReader::open_path(fixture(name)).is_err());
    }
}

#[test]
#[ignore = "requires downloaded public fixtures, run by CI"]
fn instrument_terms_use_psi_ms_accessions_and_canonical_names() {
    let path = fixture("BSA-FT-HCD.raw");
    let mut raw = RawFileReader::open_path(&path).unwrap();
    for (model, accession, canonical) in [
        ("Orbitrap Astral", "MS:1003378", "Orbitrap Astral"),
        ("Orbitrap Ascend", "MS:1003356", "Orbitrap Ascend"),
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
        ("LTQ Orbitrap Velos Pro", "MS:1003096", "Orbitrap Velos Pro"),
        ("LTQ Orbitrap Elite", "MS:1001910", "Orbitrap Elite"),
        ("LTQ Velos Pro", "MS:1003495", "Velos Pro"),
        ("TSQ Altis", "MS:1002874", "TSQ Altis"),
        ("TSQ Quantiva", "MS:1002418", "TSQ Quantiva"),
        ("TSQ Endura", "MS:1002419", "TSQ Endura"),
        ("LCQ Classic", "MS:1000168", "LCQ Classic"),
        ("LCQ Deca", "MS:1000554", "LCQ Deca"),
        ("LCQ Advantage", "MS:1000167", "LCQ Advantage"),
    ] {
        raw.instrument_model = Some(model);
        let mut source = File::open(&path).unwrap();
        let adapter = OpenTfRawSource::new(&raw, &mut source, "instrument.raw", false);
        let term = adapter.run_metadata().instrument;
        assert_eq!(term.accession, accession);
        assert_eq!(term.name, canonical);
    }
}
