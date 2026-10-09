# Local performance baseline

The rc.2 baseline uses the `reader_benchmark` release example, Rust 1.97.1 on Linux x86_64. The local acquisition is 35,488,261 bytes with 5,445 scans. No acquisition data or path is included. Each workload runs three repetitions with deterministic shuffle seed 42. Timings below are totals across three repetitions.

| Workload | Peaks per pass | Open seconds | Sequential seconds | Random seconds | Peak RSS KiB |
| --- | ---: | ---: | ---: | ---: | ---: |
| Centroid | 470,235 | 0.307286 | 0.187090 | 0.295003 | 25,984 |
| Stored profile | 2,768,240 | 0.643855 | 1.139916 | 1.031127 | 26,496 |

Peak RSS was measured externally with `/usr/bin/time`. Other builds and qualification work ran concurrently, and filesystem cache state was uncontrolled. These numbers establish a local reference, not a cross-machine guarantee or a controlled before/after speed comparison. Repeated open warms caches before scan timing. Random order measures seek patterns on the same file, not independent cold reads.

Reproduce with `cargo run --release --locked -p sage-plus-raw --example reader_benchmark -- PATH 3 42`, adding `--profile` for stored profile bins. Supply an authorized local input. The executable stops on decoding errors and verifies equal peak counts between sequential and shuffled reads.
