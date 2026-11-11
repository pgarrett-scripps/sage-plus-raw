//! Two handles on the same path against one owned handle, the pattern a
//! downstream reader uses. Run with PATH [REPEATS] [--profile].
use sage_plus_raw::{iter_spectra_from_path, try_iter_spectra, RawFileReader};
use std::{fs::File, hint::black_box, io::BufReader, time::Instant};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let path = args
        .first()
        .ok_or("usage: single_handle PATH [REPEATS] [--profile]")?;
    if args.len() > 3 || args.get(2).is_some_and(|arg| arg != "--profile") {
        return Err("unexpected argument, optional third argument is --profile".into());
    }
    let repeats: u32 = args.get(1).map(|v| v.parse()).transpose()?.unwrap_or(3);
    if repeats == 0 {
        return Err("REPEATS must be positive".into());
    }
    let profile = args.get(2).is_some_and(|v| v == "--profile");
    let scans = RawFileReader::open_path(path)?.num_scans;
    let (mut two_seconds, mut one_seconds) = (0.0, 0.0);
    let (mut two_peaks, mut one_peaks) = (0u64, 0u64);
    // One pass of each per repeat, alternating which goes first, so neither
    // path gets all of the warm page cache. Each pass includes its own opens,
    // which is what the two paths differ in.
    for repeat in 0..repeats {
        for turn in 0..2 {
            let one_handle = (turn == 0) == (repeat % 2 == 1);
            let start = Instant::now();
            if one_handle {
                for spectrum in iter_spectra_from_path(path, profile)? {
                    let spectrum = spectrum?;
                    one_peaks += spectrum.mz.len() as u64;
                    black_box(spectrum);
                }
                one_seconds += start.elapsed().as_secs_f64();
            } else {
                let raw = RawFileReader::open_path(path)?;
                let mut source = BufReader::new(File::open(path)?);
                for spectrum in try_iter_spectra(&raw, &mut source, profile) {
                    let spectrum = spectrum?;
                    two_peaks += spectrum.mz.len() as u64;
                    black_box(spectrum);
                }
                two_seconds += start.elapsed().as_secs_f64();
            }
        }
    }
    assert_eq!(two_peaks, one_peaks, "the two paths decoded differently");
    println!(
        "{{\"version\":\"{}\",\"scans\":{scans},\"repeats\":{repeats},\"profile\":{profile},\
         \"two_handle_seconds\":{two_seconds},\"one_handle_seconds\":{one_seconds},\
         \"peaks_per_pass\":{}}}",
        env!("CARGO_PKG_VERSION"),
        one_peaks / u64::from(repeats)
    );
    Ok(())
}
