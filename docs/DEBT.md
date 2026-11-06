# Technical Debt & Ratchet Ledger

Every unavoidable lint suppression, duplicate-dependency skip, echoes/coverage ratchet, or mutation exclusion in this repository is recorded here with its reason, owner, and removal condition.

## Active Exceptions

| Location | Gate / Lint | Reason | Owner | Removal Condition |
|---|---|---|---|---|
| `Cargo.toml` (`[lints.clippy]`, `[lints.rust]`, `[lints.rustdoc]`) | 17 ratcheted lints set to `"allow"` with `# N findings; when fixed: ...` comments (plus `redundant_pub_crate = "allow"`) | Legacy findings baseline captured at `rusti-graben` adoption (2026-10-30); guarded by `./mk ratchet` so no new `"allow"` entries can be added. | `@Adafede` | Fix findings per lint and flip to `"warn"` / `"deny"` or delete the line as instructed by each comment. |
| `.cargo/mutants.toml` (`exclude_re`) | 3 line:column-pinned regex exclusions in `src/reverse.rs` (`oxygen_bridge`) | Equivalent mutants in `oxygen_bridge`: a bridged ring prints only its start (`;5Ep`), and the bridge is always visited from its lower carbon first. | `@Adafede` | Permanent while `oxygen_bridge` uses two-endpoint scan; pinned to exact `line:col` so any line shift re-triggers the gate. |
| `deny.toml` (`[bans] skip`) | `skip = []` (0 skipped duplicate crates) | Zero duplicate crate versions in the dependency graph. | `@Adafede` | Keep at 0. |

## Ratchet Floors (`mk.local` & `Cargo.toml`)

- `ECHOES_MAX_APPROXIMATE`: unset in `mk.local` (uses `./mk` default floor of `0`)
- `ECHOES_MAX_REUSE`: unset in `mk.local` (uses `./mk` default floor of `0`)
- `COVERAGE_FAIL_UNDER_LINES`: unset in `mk.local` (uses `./mk` default floor of `80`)
- `Cargo.toml` lint ratchet count: `17` lints with `# N findings; when fixed:` comments
- `unsafe_code = "forbid"` (`0` `unsafe` blocks)
