//! Reproducible local benchmark. Run with PATH [REPEATS] [SEED] [--profile].
use opentfraw::{try_extract_spectrum_arrays, RawFileReader};
use std::{fs::File, hint::black_box, io::BufReader, time::Instant};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let path = args
        .first()
        .ok_or("usage: reader_benchmark PATH [REPEATS] [SEED] [--profile]")?;
    if args.len() > 4 || args.get(3).is_some_and(|arg| arg != "--profile") {
        return Err("unexpected argument, optional fourth argument is --profile".into());
    }
    let repeats: u32 = args.get(1).map(|v| v.parse()).transpose()?.unwrap_or(3);
    if repeats == 0 {
        return Err("REPEATS must be positive".into());
    }
    let seed: u64 = args.get(2).map(|v| v.parse()).transpose()?.unwrap_or(42);
    let profile = args.get(3).is_some_and(|v| v == "--profile");
    let bytes = std::fs::metadata(path)?.len();
    let start = Instant::now();
    for _ in 0..repeats {
        black_box(RawFileReader::open_path(path)?);
    }
    let open_seconds = start.elapsed().as_secs_f64();
    let raw = RawFileReader::open_path(path)?;
    let mut source = BufReader::new(File::open(path)?);
    let sequential: Vec<_> = (0..raw.num_scans).collect();
    let mut random = sequential.clone();
    let mut state = seed;
    for i in (1..random.len()).rev() {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        random.swap(i, (state % (i as u64 + 1)) as usize);
    }
    let mut run = |order: &[u32]| -> opentfraw::Result<(f64, u64)> {
        let start = Instant::now();
        let mut peaks = 0;
        for _ in 0..repeats {
            for &index in order {
                let arrays = try_extract_spectrum_arrays(&raw, &mut source, index, profile)?;
                peaks += arrays.mz.len() as u64;
                black_box((&arrays.mz, &arrays.intensity));
            }
        }
        Ok((start.elapsed().as_secs_f64(), peaks))
    };
    let (sequential_seconds, sequential_peaks) = run(&sequential)?;
    let (random_seconds, random_peaks) = run(&random)?;
    assert_eq!(sequential_peaks, random_peaks);
    println!("{{\"version\":\"{}\",\"input_bytes\":{bytes},\"scans\":{},\"repeats\":{repeats},\"seed\":{seed},\"profile\":{profile},\"open_seconds\":{open_seconds},\"sequential_seconds\":{sequential_seconds},\"random_seconds\":{random_seconds},\"peaks_per_pass\":{}}}", env!("CARGO_PKG_VERSION"), raw.num_scans, sequential_peaks / u64::from(repeats));
    Ok(())
}
