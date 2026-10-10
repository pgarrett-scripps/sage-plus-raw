# Security policy

This policy applies to the independent `sage-plus-raw` fork maintained in [pgarrett-scripps/sage-plus-raw](https://github.com/pgarrett-scripps/sage-plus-raw). It does not promise support from the upstream OpenTFRaw maintainers.

## Reporting

Report vulnerabilities privately through [this repository's security advisories](https://github.com/pgarrett-scripps/sage-plus-raw/security/advisories/new). Do not disclose exploit details in public issues, and do not send private acquisition data without permission.

Include the crate version or commit, operating system, Rust version, expected behavior and a minimal reproduction. Prefer synthetic byte sequences over sensitive RAW files. Maintainers will coordinate disclosure and attribution with the reporter. This is a volunteer-maintained project, so no acknowledgment or remediation deadline is guaranteed.

## Supported versions

Security fixes are made on the latest published release of `sage-plus-raw` on crates.io and on `main`. Earlier releases and upstream OpenTFRaw releases are outside this fork's support commitment.

## Scope

Relevant findings include malformed-input panics, excessive allocation, infinite loops, unintended file writes and reachable dependency vulnerabilities. Scientific decoding errors also matter, even when they do not create a security vulnerability.

The parser forbids unsafe Rust in its own crate. That restriction is not a guarantee that dependencies or platform libraries contain no unsafe code. Treat untrusted files as potentially hostile and apply application-level resource limits when appropriate.

## Distribution

The crate is published to crates.io only by the release workflow in `.github/workflows/release.yml`, triggered by a version tag. See `RELEASING.md`. Dependencies are checked with `cargo deny` against the RustSec advisory database on every CI run.
