//! Scratch identity check (not committed): FNV-1a over every decoded m/z and
//! intensity bit pattern, plus peak and dropped-peak counts.
//! Run with PATH [SCANS] [--profile].
use sage_plus_raw::{try_extract_spectrum_into, RawFileReader};
use std::{fs::File, io::BufReader};

fn fold(hash: &mut u64, bytes: &[u8]) {
    for byte in bytes {
        *hash ^= u64::from(*byte);
        *hash = hash.wrapping_mul(0x100_0000_01b3);
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let path = args.first().ok_or("usage: zz_checksum PATH [SCANS] [--profile]")?;
    let raw = RawFileReader::open_path(path)?;
    let scans: u32 = args
        .get(1)
        .map(|v| v.parse())
        .transpose()?
        .unwrap_or(raw.num_scans)
        .min(raw.num_scans);
    let profile = args.iter().any(|v| v == "--profile");
    let mut source = BufReader::new(File::open(path)?);
    let mut mz = Vec::new();
    let mut intensity = Vec::new();
    let mut hash = 0xcbf2_9ce4_8422_2325u64;
    let (mut peaks, mut dropped) = (0u64, 0u64);
    for index in 0..scans {
        let header =
            try_extract_spectrum_into(&raw, &mut source, index, profile, &mut mz, &mut intensity)?;
        peaks += mz.len() as u64;
        dropped += header.dropped_peaks as u64;
        fold(&mut hash, &(mz.len() as u64).to_le_bytes());
        for value in &mz {
            fold(&mut hash, &value.to_bits().to_le_bytes());
        }
        for value in &intensity {
            fold(&mut hash, &value.to_bits().to_le_bytes());
        }
        fold(&mut hash, &(header.dropped_peaks as u64).to_le_bytes());
        fold(&mut hash, &format!("{:?}", header.scan_mode).into_bytes());
    }
    println!(
        "{{\"scans\":{scans},\"profile\":{profile},\"peaks\":{peaks},\"dropped\":{dropped},\"hash\":\"{hash:016x}\"}}"
    );
    Ok(())
}
