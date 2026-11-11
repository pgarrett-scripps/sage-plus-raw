//! Fresh-allocation against reused-buffer decoding, the pattern a downstream
//! reader loop uses. Run with PATH [REPEATS] [SCANS] [--profile].
use sage_plus_raw::{try_extract_spectrum, try_extract_spectrum_into, RawFileReader};
use std::{fs::File, hint::black_box, io::BufReader, time::Instant};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let path = args
        .first()
        .ok_or("usage: reuse_buffers PATH [REPEATS] [SCANS] [--profile]")?;
    if args.len() > 4 || args.get(3).is_some_and(|arg| arg != "--profile") {
        return Err("unexpected argument, optional fourth argument is --profile".into());
    }
    let repeats: u32 = args.get(1).map(|v| v.parse()).transpose()?.unwrap_or(3);
    if repeats == 0 {
        return Err("REPEATS must be positive".into());
    }
    let raw = RawFileReader::open_path(path)?;
    let scans: u32 = args
        .get(2)
        .map(|v| v.parse())
        .transpose()?
        .unwrap_or(raw.num_scans)
        .min(raw.num_scans);
    let profile = args.get(3).is_some_and(|v| v == "--profile");
    let mut source = BufReader::new(File::open(path)?);
    let mut mz = Vec::new();
    let mut intensity = Vec::new();
    let (mut fresh_seconds, mut reused_seconds) = (0.0, 0.0);
    let (mut fresh_peaks, mut reused_peaks) = (0u64, 0u64);
    // One pass of each per repeat, alternating which goes first, so neither
    // path gets all of the warm page cache.
    for repeat in 0..repeats {
        for turn in 0..2 {
            let reused = (turn == 0) == (repeat % 2 == 1);
            let start = Instant::now();
            if reused {
                for index in 0..scans {
                    let header = try_extract_spectrum_into(
                        &raw,
                        &mut source,
                        index,
                        profile,
                        &mut mz,
                        &mut intensity,
                    )?;
                    reused_peaks += mz.len() as u64;
                    black_box((&header, &mz, &intensity));
                }
                reused_seconds += start.elapsed().as_secs_f64();
            } else {
                for index in 0..scans {
                    let record = try_extract_spectrum(&raw, &mut source, index, profile)?;
                    fresh_peaks += record.mz.len() as u64;
                    black_box(record);
                }
                fresh_seconds += start.elapsed().as_secs_f64();
            }
        }
    }
    assert_eq!(
        fresh_peaks, reused_peaks,
        "the two paths decoded differently"
    );
    println!(
        "{{\"version\":\"{}\",\"scans\":{scans},\"repeats\":{repeats},\"profile\":{profile},\
         \"fresh_seconds\":{fresh_seconds},\"reused_seconds\":{reused_seconds},\
         \"peaks_per_pass\":{}}}",
        env!("CARGO_PKG_VERSION"),
        fresh_peaks / u64::from(repeats)
    );
    Ok(())
}
