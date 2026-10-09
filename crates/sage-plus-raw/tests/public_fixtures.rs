//! Independent ProteoWizard acquisition regressions, required by CI.
use sage_plus_raw::{Activation, RawFileReader, ScanMode};
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
        let meta = sage_plus_raw::scan_metadata(raw, 0).unwrap();
        assert_eq!(
            meta.precursor.as_ref().unwrap().activation,
            Some(activation)
        );
        assert!(meta.precursor.unwrap().master_scan_number.is_none());
        assert!(meta.filter.unwrap().contains(token));
    }
}

#[test]
#[ignore = "requires downloaded public fixtures, run by CI"]
fn strict_paths_propagate_decode_errors_and_label_centroids() {
    let path = fixture("BSA-FT-HCD.raw");
    let raw = RawFileReader::open_path(&path).unwrap();
    assert_eq!(
        sage_plus_raw::scan_metadata(&raw, 0).unwrap().scan_mode,
        Some(ScanMode::Profile)
    );
    let spectrum =
        sage_plus_raw::try_extract_spectrum(&raw, &mut File::open(&path).unwrap(), 0, false)
            .unwrap();
    assert_eq!(spectrum.scan_mode, Some(ScanMode::Centroid));
    assert!(!spectrum.mz.is_empty());
    assert!(sage_plus_raw::try_extract_spectrum(
        &raw,
        &mut Cursor::new(Vec::<u8>::new()),
        0,
        false
    )
    .is_err());
    assert!(
        sage_plus_raw::try_extract_spectrum(&raw, &mut File::open(&path).unwrap(), 1, false)
            .is_err()
    );
    let results: Vec<_> =
        sage_plus_raw::try_iter_spectra(&raw, &mut Cursor::new(Vec::<u8>::new()), false).collect();
    assert_eq!(results.len(), raw.num_scans as usize);
    assert!(results.iter().all(Result::is_err));
}

#[test]
#[ignore = "requires downloaded public fixtures, run by CI"]
fn unknown_activation_tags_stay_unknown() {
    let path = fixture("BSA-FT-HCD.raw");
    let mut raw = RawFileReader::open_path(&path).unwrap();
    raw.scan_events[0].reactions[0].unknown_long1 = 999;
    assert!(sage_plus_raw::scan_metadata(&raw, 0)
        .unwrap()
        .precursor
        .unwrap()
        .activation
        .is_none());
    assert!(!raw.scan_filter(1).unwrap().contains('@'));
    raw.scan_events[0].reactions[0].unknown_long1 = 11;
    raw.scan_events[0].preamble.bytes[6] = 3;
    let meta = sage_plus_raw::scan_metadata(&raw, 0).unwrap();
    assert_eq!(meta.ms_level, 3);
    assert!(meta.precursor.is_some());
}

#[test]
#[ignore = "requires downloaded public fixtures, run by CI"]
fn profile_requires_calibration_and_strict_metadata_requires_an_event() {
    let path = fixture("BSA-FT-HCD.raw");
    let mut raw = RawFileReader::open_path(&path).unwrap();
    let profile =
        sage_plus_raw::try_extract_spectrum(&raw, &mut File::open(&path).unwrap(), 0, true)
            .unwrap();
    assert_eq!(profile.scan_mode, Some(ScanMode::Profile));
    assert!(!profile.mz.is_empty());
    raw.scan_events[0].coefficients.clear();
    assert!(
        sage_plus_raw::try_extract_spectrum(&raw, &mut File::open(&path).unwrap(), 0, true)
            .is_err()
    );
    raw.scan_events.clear();
    assert!(
        sage_plus_raw::try_extract_spectrum(&raw, &mut File::open(&path).unwrap(), 0, false)
            .is_err()
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
    // The precursor describes the first multiplexed reaction.
    let precursor = sage_plus_raw::scan_metadata(&raw, 0)
        .unwrap()
        .precursor
        .unwrap();
    assert!((precursor.target_mz.unwrap() - expected[0]).abs() < 1e-9);
    let raw =
        RawFileReader::open_path(fixture("IsolationMzOffset-ReportedMassOffset.raw")).unwrap();
    let precursor = sage_plus_raw::scan_metadata(&raw, 0)
        .unwrap()
        .precursor
        .unwrap();
    assert!((precursor.target_mz.unwrap() - 379.7113342285156).abs() < 1e-9);
    assert!((precursor.selected_mz.unwrap() - 384.7153374006682).abs() < 1e-9);
    // Variable-length event streams open and decode every event in step.
    for (name, levels) in [
        ("IT-HCD-SPS.raw", &[3u8][..]),
        ("source_cid_test_3scans.raw", &[2, 2, 2][..]),
    ] {
        let raw = RawFileReader::open_path(fixture(name)).unwrap();
        let decoded: Vec<u8> = raw
            .scan_events
            .iter()
            .map(|event| event.preamble.bytes[6])
            .collect();
        assert_eq!(decoded, levels, "{name}");
    }
}

#[test]
#[ignore = "requires normalized public references, run by CI"]
fn full_precision_arrays_match_independent_proteowizard_references() {
    use std::io::Read;
    fn u32_value(reader: &mut impl Read) -> u32 {
        let mut bytes = [0; 4];
        reader.read_exact(&mut bytes).unwrap();
        u32::from_le_bytes(bytes)
    }
    fn f64_value(reader: &mut impl Read) -> f64 {
        let mut bytes = [0; 8];
        reader.read_exact(&mut bytes).unwrap();
        f64::from_le_bytes(bytes)
    }
    for name in [
        "BSA-FT-ETD",
        "BSA-FT-HCD",
        "FT-HCD-MSX",
        "IsolationMzOffset-ReportedMassOffset",
    ] {
        let path = fixture(&format!("{name}.raw"));
        let raw = RawFileReader::open_path(&path).unwrap();
        for profile in [false, true] {
            let suffix = if profile { "" } else { "-centroid" };
            let mut reference = File::open(fixture(&format!("{name}{suffix}.arrays"))).unwrap();
            assert_eq!(u32_value(&mut reference), raw.num_scans, "{name}");
            for _ in 0..raw.num_scans {
                let idx = u32_value(&mut reference);
                let level = u32_value(&mut reference);
                let rt = f64_value(&mut reference);
                let count = u32_value(&mut reference) as usize;
                let arrays = sage_plus_raw::try_extract_spectrum_arrays(
                    &raw,
                    &mut File::open(&path).unwrap(),
                    idx,
                    profile,
                )
                .unwrap();
                assert_eq!(arrays.metadata.ms_level, level, "{name} scan {idx}");
                assert!(
                    (arrays.metadata.retention_time_min - rt).abs() < 1e-6,
                    "{name} scan {idx} RT"
                );
                let pairs: Vec<_> = (0..count)
                    .map(|_| (f64_value(&mut reference), f64_value(&mut reference)))
                    .filter(|(_, intensity)| !profile || *intensity != 0.0)
                    .collect();
                assert_eq!(
                    arrays.mz.len(),
                    pairs.len(),
                    "{name} scan {idx} profile={profile}"
                );
                for (peak, (mz, intensity)) in pairs.into_iter().enumerate() {
                    assert!(
                        (arrays.mz[peak] - mz).abs() <= f64::from(f32::EPSILON) * mz.abs(),
                        "{name} scan {idx} peak {peak} profile={profile} mz {} != {mz}",
                        arrays.mz[peak]
                    );
                    assert!((arrays.intensity[peak] - intensity).abs() <= f64::from(f32::EPSILON) * intensity.abs().max(1.0), "{name} scan {idx} peak {peak} profile={profile} intensity {} != {intensity}", arrays.intensity[peak]);
                }
            }
        }
    }
}

#[test]
#[ignore = "requires downloaded public fixtures, run by CI"]
fn failures_retain_scan_offsets_categories_and_io_sources() {
    use sage_plus_raw::{ErrorKind, Operation};
    let path = fixture("BSA-FT-HCD.raw");
    let raw = RawFileReader::open_path(&path).unwrap();
    let error = sage_plus_raw::try_extract_spectrum_arrays(&raw, &mut Cursor::new([]), 0, false)
        .unwrap_err();
    assert_eq!(error.kind(), ErrorKind::Corrupt);
    let context = error.context().unwrap();
    assert_eq!(context.operation, Operation::DecodeSpectrum);
    assert_eq!(context.scan_index, Some(0));
    assert_eq!(context.scan_number, Some(1));
    assert!(context.byte_offset.is_some_and(|offset| offset > 0));
    assert!(matches!(
        error.root_cause(),
        sage_plus_raw::Error::UnexpectedEof { .. } | sage_plus_raw::Error::Io(_)
    ));
    let invalid = sage_plus_raw::try_extract_spectrum_arrays(
        &raw,
        &mut File::open(&path).unwrap(),
        u32::MAX,
        false,
    )
    .unwrap_err();
    assert_eq!(invalid.kind(), ErrorKind::InvalidInput);
    assert_eq!(invalid.context().unwrap().scan_index, Some(u32::MAX));
}

#[test]
#[ignore = "requires downloaded public fixtures, run by CI"]
fn faims_requires_explicit_enablement_and_preserves_active_zero() {
    use sage_plus_raw::generic_data::GenericValue;
    let mut raw = RawFileReader::open_path(fixture("BSA-FT-HCD.raw")).unwrap();
    for (enabled, voltage, expected) in [
        (Some(GenericValue::Bool(false)), 0.0, None),
        (Some(GenericValue::Bool(false)), -45.0, None),
        (Some(GenericValue::Bool(true)), 0.0, Some(0.0)),
        (Some(GenericValue::Bool(true)), -45.0, Some(-45.0)),
        (Some(GenericValue::Bool(true)), f64::NAN, None),
        (Some(GenericValue::String(" On ".into())), 0.0, Some(0.0)),
        (Some(GenericValue::String("Off".into())), -45.0, None),
        (Some(GenericValue::String("not on".into())), -45.0, None),
        (None, -45.0, None),
    ] {
        let params = &mut raw.scan_parameters[0].values;
        params.retain(|(key, _)| key != "FAIMS CV:" && key != "FAIMS Voltage On:");
        params.push(("FAIMS CV:".into(), GenericValue::Float64(voltage)));
        if let Some(enabled) = enabled {
            params.push(("FAIMS Voltage On:".into(), enabled));
        }
        assert_eq!(
            sage_plus_raw::scan_metadata(&raw, 0).unwrap().faims_cv,
            expected
        );
    }
}
