# Security policy

This policy applies to the independent `opentfraw-reader` fork maintained in [pgarrett-scripps/opentfraw](https://github.com/pgarrett-scripps/opentfraw). It does not promise support from the upstream OpenTFRaw maintainers.

## Reporting

While the repository is private, collaborators should contact its owner through their established private project channel. Do not disclose exploit details in public issues or send private acquisition data without permission.

Before making the repository public, maintainers must enable and verify GitHub private vulnerability reporting. The intended public reporting location is [this repository's security advisories](https://github.com/pgarrett-scripps/opentfraw/security/advisories). Its availability must be checked during the launch gate. This document does not claim the public reporting service is already enabled.

Include the commit or package version, operating system, Rust version, expected behavior and a minimal reproduction. Prefer synthetic byte sequences over sensitive RAW files. Maintainers will coordinate disclosure and attribution with the reporter. No acknowledgment or remediation deadline is guaranteed during release preparation.

## Scope and supported versions

The current launch candidate and the most recent private release are evaluated for security fixes. Earlier upstream releases are outside this fork's support commitment. No public crate release currently exists under the new package identity.

Relevant findings include malformed-input panics, excessive allocation, infinite loops, unintended file writes and reachable dependency vulnerabilities. Scientific decoding errors also matter, even when they do not create a security vulnerability.

The parser forbids unsafe Rust in its own crate. That restriction is not a guarantee that dependencies or platform libraries contain no unsafe code. Treat untrusted files as potentially hostile and apply application-level resource limits when appropriate.

## Distribution

The real Rust manifests retain `publish = false`. Python packaging retains a private distribution classifier. CI runs registry dry-runs only inside disposable copies. There are no registry upload, website deployment or visibility-changing workflows.
