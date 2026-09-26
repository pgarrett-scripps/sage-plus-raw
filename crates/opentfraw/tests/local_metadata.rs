use opentfraw::RawFileReader;
use std::path::PathBuf;

#[test]
#[ignore = "requires THERMO_LARGE_RAW and THERMO_PUBLIC_CORPUS local fixture paths"]
fn local_metadata_regressions() {
    let large = std::env::var("THERMO_LARGE_RAW").unwrap();
    let raw = RawFileReader::open_path(&large).unwrap();
    assert_eq!(raw.num_scans, 5445);
    let mut source = std::fs::File::open(&large).unwrap();
    let mut scans = 0;
    let mut peaks = 0;
    for spectrum in opentfraw::try_iter_spectra(&raw, &mut source, false) {
        let spectrum = spectrum.unwrap();
        scans += 1;
        peaks += spectrum.mz.len();
        assert_eq!(spectrum.scan_mode, Some(opentfraw::ScanMode::Centroid));
        assert_eq!(spectrum.faims_cv, None);
        assert_eq!(
            raw.scan_params(spectrum.scan_number)
                .unwrap()
                .faims_voltage_on(),
            Some(false)
        );
    }
    assert_eq!(scans, 5445);
    assert_eq!(peaks, 470235);
    assert_eq!(
        raw.scan_events
            .iter()
            .filter(|e| e.preamble.bytes[6] == 2)
            .count(),
        2000
    );
    assert!(raw.scan_events.iter().all(|e| e.coefficients.len() == 7));
    let corpus = PathBuf::from(std::env::var("THERMO_PUBLIC_CORPUS").unwrap());
    let path = corpus.join("FT-HCD-MSX.raw");
    let raw = RawFileReader::open_path(&path).unwrap();
    let expected = [522.987607, 698.067188, 502.978512, 768.099021];
    assert_eq!(raw.scan_events[0].reactions.len(), expected.len());
    for (r, mz) in raw.scan_events[0].reactions.iter().zip(expected) {
        assert!((r.precursor_mz - mz).abs() < 1e-9);
        assert_eq!(r.unknown_double, 5.0);
    }
    assert!(opentfraw::scan_metadata(&raw, 0)
        .unwrap()
        .precursor
        .is_none());
    let mut output = Vec::new();
    assert!(opentfraw::write_mzml(
        &raw,
        &mut std::fs::File::open(path).unwrap(),
        &mut output,
        "msx.raw",
        false
    )
    .is_err());
    assert!(output.is_empty());
    assert!(RawFileReader::open_path(corpus.join("source_cid_test_3scans.raw")).is_err());
}
