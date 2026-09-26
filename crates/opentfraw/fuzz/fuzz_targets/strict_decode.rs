#![no_main]
use libfuzzer_sys::fuzz_target;
use opentfraw::{try_extract_spectrum_arrays, RawFileReader};
use std::io::Cursor;

fuzz_target!(|data: &[u8]| {
    if let Ok(raw) = RawFileReader::open(Cursor::new(data)) {
        let mut source = Cursor::new(data);
        for index in 0..raw.num_scans.min(4) {
            for profile in [false, true] {
                let _ = try_extract_spectrum_arrays(&raw, &mut source, index, profile);
            }
        }
    }
});
