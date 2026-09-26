//! Print event metadata for qualification against independent acquisition references.
use opentfraw::RawFileReader;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let path = std::env::args()
        .nth(1)
        .ok_or("usage: review_metadata file.raw")?;
    let raw = RawFileReader::open_path(&path)?;
    println!(
        "version={} trailer={} params={}",
        raw.header.version, raw.run_header.scan_trailer_addr, raw.run_header.scan_params_addr
    );
    let mut levels = std::collections::BTreeMap::new();
    for (i, event) in raw.scan_events.iter().enumerate() {
        let level = event.preamble.bytes.get(6).copied().unwrap_or(0);
        *levels.entry(level).or_insert(0usize) += 1;
        if i < 3 {
            println!(
                "event={} level={} reactions={:?} coefficients={:?}",
                i + 1,
                level,
                event.reactions,
                event.coefficients
            );
        }
    }
    println!(
        "scans={} events={} levels={levels:?}",
        raw.num_scans,
        raw.scan_events.len()
    );
    Ok(())
}
