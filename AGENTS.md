# AGENTS.md

This repository uses the [rusti-graben](https://github.com/Adafede/rusti-graben) quality gates.

- **Single entrypoint**: run `./mk <task>` (`./mk --list` shows them). CI calls the same tasks, so a green `./mk ci` locally means a green CI.
- **Never weaken a gate.** Record unavoidable exceptions where they are reviewed:
  `#[expect(lint, reason = "...")]` (not `#[allow]`), `// echoes: allow` above an intentional copy,
  a commented `.cargo/mutants.toml` exclusion, a `deny.toml` `skip` entry with a reason.
- **Pinned tools**: versions live in `mk.versions`; `./mk setup` installs them. Do not install other versions.
- **Commits** follow Conventional Commits (`feat:`, `fix:`, `chore(deps):`, ...); `./mk commits` checks the branch.
- **Before finishing a change**: `./mk ci`. New code needs tests that kill its mutants (`./mk mutants-diff`).
