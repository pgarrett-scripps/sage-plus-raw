//! Stream every spectrum and propagate the first decoding error.
use sage_plus_raw::{try_iter_spectra, RawFileReader};
use std::{fs::File, io::BufReader};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let path = std::env::args()
        .nth(1)
        .ok_or("usage: cargo run --example strict_read -- acquisition.raw")?;
    let raw = RawFileReader::open_path(&path)?;
    let mut source = BufReader::new(File::open(&path)?);
    for spectrum in try_iter_spectra(&raw, &mut source, false) {
        let spectrum = spectrum?;
        println!(
            "scan={} ms_level={} peaks={}",
            spectrum.scan_number,
            spectrum.ms_level,
            spectrum.mz.len()
        );
    }
    Ok(())
}
