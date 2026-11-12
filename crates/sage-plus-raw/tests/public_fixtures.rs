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
        // Decoding reads the compact trailer row, taken out of the record
        // while the file was opened, so an edited record reaches it the way
        // the open loop puts it there.
        raw.scan_trailers[0] = sage_plus_raw::ScanTrailer::from_params(&sage_plus_raw::ScanParams(
            &raw.scan_parameters[0],
        ));
        assert_eq!(
            sage_plus_raw::scan_metadata(&raw, 0).unwrap().faims_cv,
            expected
        );
    }
}

/// Every acquisition in the public corpus, decoded twice: once per scan into
/// a fresh record, once into two buffers reused across the whole file. The
/// reused path must produce the identical record, filter string and error.
#[test]
#[ignore = "requires downloaded public fixtures, run by CI"]
fn scan_index_entries_are_eighty_eight_bytes_apart_in_version_66_files() {
    // The scan index is read back to back from `scan_index_addr`, so the
    // entry size is the stride. Every public fixture is version 66, and in
    // each the bytes between the scan index and the trailer stream that
    // follows it are exactly `num_scans * 88`. 92 bytes, which this struct's
    // documentation claimed before, would overrun the trailer.
    for name in [
        "BSA-FT-ETD.raw",
        "BSA-FT-HCD.raw",
        "FT-HCD-MSX.raw",
        "IsolationMzOffset-ReportedMassOffset.raw",
        "IT-HCD-SPS.raw",
        "source_cid_test_3scans.raw",
    ] {
        let raw = RawFileReader::open_path(fixture(name)).unwrap();
        assert_eq!(raw.version, 66, "{name}");
        let header = &raw.run_header;
        let span = header
            .scan_trailer_addr
            .checked_sub(header.scan_index_addr)
            .unwrap_or_else(|| panic!("{name}: trailer does not follow the scan index"));
        assert_eq!(span, u64::from(raw.num_scans) * 88, "{name}");
        // The entries decoded at that stride are coherent: sequential
        // indices, non-decreasing start times, and each scan's packet ending
        // where the next one begins.
        for (position, entry) in raw.scan_index.iter().enumerate() {
            assert_eq!(u64::from(entry.index), position as u64, "{name}");
        }
        for pair in raw.scan_index.windows(2) {
            assert!(pair[1].start_time >= pair[0].start_time, "{name}");
            assert_eq!(
                pair[1].offset,
                pair[0].offset + u64::from(pair[0].data_size),
                "{name}"
            );
        }
    }
}

#[test]
#[ignore = "requires downloaded public fixtures, run by CI"]
fn reused_buffers_decode_what_fresh_allocations_do() {
    for name in [
        "BSA-FT-ETD.raw",
        "BSA-FT-HCD.raw",
        "FT-HCD-MSX.raw",
        "IT-HCD-SPS.raw",
        "IsolationMzOffset-ReportedMassOffset.raw",
        "source_cid_test_3scans.raw",
    ] {
        let path = fixture(name);
        let raw = RawFileReader::open_path(&path).unwrap();
        let mut mz = Vec::new();
        let mut intensity = Vec::new();
        let mut filter = String::new();
        for include_profile in [false, true] {
            let mut fresh = File::open(&path).unwrap();
            let mut reused = File::open(&path).unwrap();
            for idx in 0..raw.num_scans {
                let expected =
                    sage_plus_raw::try_extract_spectrum(&raw, &mut fresh, idx, include_profile);
                let header = sage_plus_raw::try_extract_spectrum_into(
                    &raw,
                    &mut reused,
                    idx,
                    include_profile,
                    &mut mz,
                    &mut intensity,
                );
                match (expected, header) {
                    (Ok(record), Ok(header)) => {
                        assert_eq!(record.mz, mz, "{name} scan index {idx} m/z");
                        assert_eq!(record.intensity, intensity, "{name} scan index {idx}");
                        let rebuilt = header.into_record(mz.clone(), intensity.clone());
                        assert_eq!(
                            format!("{record:?}"),
                            format!("{rebuilt:?}"),
                            "{name} scan index {idx} record"
                        );
                    }
                    (Err(expected), Err(found)) => {
                        assert_eq!(expected.to_string(), found.to_string());
                        assert!(mz.is_empty() && intensity.is_empty());
                    }
                    (expected, found) => {
                        panic!("{name} scan index {idx}: {expected:?} vs {found:?}")
                    }
                }

                let scan_number = raw.run_header.sample_info.first_scan_number + idx;
                let rendered = raw.scan_filter_into(scan_number, &mut filter);
                match raw.scan_filter(scan_number) {
                    Some(owned) => {
                        assert!(rendered, "{name} scan {scan_number}");
                        assert_eq!(owned, filter, "{name} scan {scan_number} filter");
                    }
                    None => {
                        assert!(!rendered && filter.is_empty(), "{name} scan {scan_number}");
                    }
                }
            }
        }
    }
}

#[test]
#[ignore = "requires downloaded public fixtures, run by CI"]
fn one_handle_decodes_what_two_handles_do() {
    // The sequential reader holds the handle it parsed the file with, so a
    // caller does not open the same path twice. Same bytes in, same records
    // out, scan for scan.
    for name in [
        "BSA-FT-ETD.raw",
        "BSA-FT-HCD.raw",
        "FT-HCD-MSX.raw",
        "IT-HCD-SPS.raw",
        "IsolationMzOffset-ReportedMassOffset.raw",
        "source_cid_test_3scans.raw",
    ] {
        let path = fixture(name);
        let raw = RawFileReader::open_path(&path).unwrap();
        for include_profile in [false, true] {
            let mut two_handles = File::open(&path).unwrap();
            let expected = sage_plus_raw::try_iter_spectra(&raw, &mut two_handles, include_profile);
            let found = sage_plus_raw::iter_spectra_from_path(&path, include_profile).unwrap();
            assert_eq!(found.len(), raw.num_scans as usize, "{name}");
            let mut scans = 0;
            for (idx, (expected, found)) in expected.zip(found).enumerate() {
                scans += 1;
                match (expected, found) {
                    (Ok(expected), Ok(found)) => assert_eq!(
                        format!("{expected:?}"),
                        format!("{found:?}"),
                        "{name} scan index {idx}"
                    ),
                    (Err(expected), Err(found)) => {
                        assert_eq!(
                            expected.to_string(),
                            found.to_string(),
                            "{name} index {idx}"
                        )
                    }
                    (expected, found) => {
                        panic!("{name} scan index {idx}: {expected:?} vs {found:?}")
                    }
                }
            }
            assert_eq!(scans, raw.num_scans, "{name}");
        }

        // Reading by index through the owned handle, and into reused
        // buffers, agrees with the same scan read the existing way.
        let mut reader = sage_plus_raw::SpectrumReader::open_path(&path).unwrap();
        assert_eq!(reader.num_scans(), raw.num_scans);
        assert_eq!(reader.path(), path);
        assert_eq!(reader.raw().num_scans, raw.num_scans);
        let mut handle = File::open(&path).unwrap();
        let mut mz = Vec::new();
        let mut intensity = Vec::new();
        // Back to front: the buffer suits a forward read, the index must not.
        for idx in (0..raw.num_scans).rev() {
            let expected = sage_plus_raw::try_extract_spectrum(&raw, &mut handle, idx, false);
            let header = reader.read_spectrum_into(idx, false, &mut mz, &mut intensity);
            match (expected, header) {
                (Ok(expected), Ok(header)) => {
                    let found = header.into_record(mz.clone(), intensity.clone());
                    assert_eq!(
                        format!("{expected:?}"),
                        format!("{found:?}"),
                        "{name} {idx}"
                    );
                }
                (Err(expected), Err(found)) => {
                    assert_eq!(expected.to_string(), found.to_string(), "{name} {idx}");
                    assert!(mz.is_empty() && intensity.is_empty());
                }
                (expected, found) => panic!("{name} scan index {idx}: {expected:?} vs {found:?}"),
            }
        }
    }
}

#[test]
#[ignore = "requires downloaded public fixtures, run by CI"]
fn dropping_the_per_scan_records_decodes_the_same_spectra() {
    // `OpenOptions::spectra_only` leaves the per-scan trailer and
    // instrument-log records unparsed, which is what makes a concurrent
    // five-file search hold hundreds of megabytes less. Everything decoding
    // reads from a trailer record is taken out while the file is opened, so
    // every scan's metadata, filter string and peak arrays must come out
    // byte for byte as they do with the records kept.
    for name in [
        "BSA-FT-ETD.raw",
        "BSA-FT-HCD.raw",
        "FT-HCD-MSX.raw",
        "IsolationMzOffset-ReportedMassOffset.raw",
        "IT-HCD-SPS.raw",
        "source_cid_test_3scans.raw",
    ] {
        let path = fixture(name);
        let kept = RawFileReader::open_path(&path).unwrap();
        let dropped =
            RawFileReader::open_path_with(&path, sage_plus_raw::OpenOptions::spectra_only())
                .unwrap();
        assert!(!kept.scan_parameters.is_empty(), "{name}");
        assert!(dropped.scan_parameters.is_empty(), "{name}");
        assert!(dropped.inst_log.is_empty(), "{name}");
        // The compact table is there either way, and is the same table.
        assert_eq!(kept.scan_trailers.len(), kept.num_scans as usize, "{name}");
        assert_eq!(kept.scan_trailers, dropped.scan_trailers, "{name}");
        assert_eq!(
            dropped.scan_parameters_header.fields.len(),
            kept.scan_parameters_header.fields.len(),
            "{name}"
        );

        for scan_number in kept.run_header.sample_info.first_scan_number
            ..=kept.run_header.sample_info.last_scan_number
        {
            assert_eq!(
                kept.scan_filter(scan_number),
                dropped.scan_filter(scan_number),
                "{name} scan {scan_number}"
            );
            assert_eq!(
                kept.source_cid_energy_ev(scan_number),
                dropped.source_cid_energy_ev(scan_number),
                "{name} scan {scan_number}"
            );
            assert_eq!(
                kept.source_cid_energy_source(scan_number),
                dropped.source_cid_energy_source(scan_number),
                "{name} scan {scan_number}"
            );
        }
        for idx in 0..kept.num_scans {
            assert_eq!(
                format!("{:?}", sage_plus_raw::scan_metadata(&kept, idx)),
                format!("{:?}", sage_plus_raw::scan_metadata(&dropped, idx)),
                "{name} scan index {idx}"
            );
        }
        for include_profile in [false, true] {
            let mut with_records = sage_plus_raw::SpectrumReader::open_path(&path).unwrap();
            let mut without = sage_plus_raw::SpectrumReader::open_path_with(
                &path,
                sage_plus_raw::OpenOptions::spectra_only(),
            )
            .unwrap();
            for idx in 0..kept.num_scans {
                let expected = with_records.read_spectrum(idx, include_profile);
                let found = without.read_spectrum(idx, include_profile);
                match (expected, found) {
                    (Ok(expected), Ok(found)) => assert_eq!(
                        format!("{expected:?}"),
                        format!("{found:?}"),
                        "{name} scan index {idx}"
                    ),
                    (Err(expected), Err(found)) => assert_eq!(
                        expected.to_string(),
                        found.to_string(),
                        "{name} scan index {idx}"
                    ),
                    (expected, found) => {
                        panic!("{name} scan index {idx}: {expected:?} vs {found:?}")
                    }
                }
            }
        }
    }
}
