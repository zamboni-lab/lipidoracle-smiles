# AGENTS.md

`lipidoracle-smiles` (`lipid_notation`) is a pure-Rust library that converts Shorthand2020 lipid names to SMILES/CXSMILES and back without inventing undetermined structural details. Fully resolved lipids emit plain SMILES, whereas unresolved double-bond positions (`Sg:` + `constrain(...)`), unlocalized modifications (`m:`), and unresolved sn-assignments (`$snN$` + `swappable(...)`) are faithfully encoded in extended CXSMILES layers and verified on round-trip recognition.

## Architecture Map

| Crate / Module | Responsibility | Key Types / Entry Points |
|---|---|---|
| `src/lib.rs` (`lipid_notation`) | Public bidirectional Shorthand2020 <-> SMILES/CXSMILES API | `name2smiles`, `smiles2name`, `name2structure`, `class_needs_multi_chain`, `canonicalize`, `smiles_for_depiction`, `smiles_expand`, `LipidStructure`, `ChainAtoms` |
| `src/nomenclature.rs` | Shorthand2020 headgroup, chain, double-bond, and modification grammar | `nomenclature` |
| `src/forward.rs` | Forward compilation from parsed Shorthand2020 names to SMILES/CXSMILES and atom maps | `generate_smiles`, `generate_structure`, `smiles_expand`, `LipidStructure`, `ChainAtoms` |
| `src/reverse.rs` | Reverse recognition from canonicalized SMILES/CXSMILES back to Shorthand2020 with round-trip proof | `parse_smiles` |
| `src/cxsmiles.rs` | CXSMILES trailer parsing, canonicalization, and graph-coupled reindexing | `canonicalize` |
| `src/depiction.rs` | Deterministic 2D depiction expansion of `m:` position variation | `smiles_for_depiction` |
| `src/consensus.rs` | Bracketed consensus tail (`[DB ...]`) encoding and recovery via `dbPos`/`mPos` trailer tokens | `consensus` |

See [`docs/ARCHITECTURE.md`](./docs/ARCHITECTURE.md) for module boundaries and pipeline stages.

## Commands

Always invoke checks and builds through `./mk` so local runs match pre-push hooks and CI at pinned tool versions (`mk.versions`):

```bash
./mk setup                     # Install pinned toolchains and tools (once)
./mk check                     # Fast inner loop (<10s): fmt, tombi, typos, clippy
./mk fix                       # Auto-fix formatting, TOML, typos, and clippy
./mk test-fast                 # Unit + integration tests (nextest / cargo test, no doctests)
./mk test-fast -E 'test(name)' # Run one targeted test via nextest
./mk test                      # Full test suite including doctests and README examples
./mk watch                     # Continuous feedback via bacon / cargo-watch
./mk mutants-diff              # Mutation testing on branch diff against base
./mk ci                        # All blocking PR gates (must pass before finishing)
./mk ci-nightly                # Nightly gates: mutants, coverage (cov), powerset (hack), miri, outdated
```

Set `NO_COLOR=1` or `CI=1` for machine-readable, uncolored output.

## Definition of Done

1. `./mk ci` passes locally with zero errors.
2. Every new or modified code path has deterministic tests that kill its mutants (`./mk mutants-diff`).
3. Public items have rustdoc comments with compiling examples where applicable (`./mk doc`).
4. All commits on the branch follow Conventional Commits (`feat:`, `fix:`, `refactor:`, `docs:`, `test:`, `chore:`; verified by `./mk commits`).

## Conventions & Hard Prohibitions

- **Single entrypoint**: Always invoke `./mk <task>`. Never edit `mk`, `mk.versions`, or `.github/actions/setup/action.yml`; put repo-specific settings in `mk.local`.
- **Never weaken a gate**:
  - No `unsafe` (`unsafe_code = "forbid"`).
  - No `unwrap()`, `expect()`, `panic!()`, `todo!()`, or `dbg!()` outside tests (`unwrap_used` and `expect_used` are ratcheted and may only decrease).
  - No bare `#[allow(...)]`; use `#[expect(lint, reason = "...")]` so the attribute fails once no longer needed.
  - Record any unavoidable legacy exception in [`docs/DEBT.md`](./docs/DEBT.md) alongside its inline marker (`// echoes: allow`, `// refquiem: allow(name)`, `.cargo/mutants.toml` comment, or `deny.toml` `skip` reason).
- **Deterministic tests**: No live network calls, wall-clock dependence, or unseeded randomness in unit/integration tests.
