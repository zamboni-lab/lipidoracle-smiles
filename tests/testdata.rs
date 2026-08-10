//! Corpus tests driven by the CSVs in `testdata/`.
//!
//! Two complementary checks live here:
//!
//! * **Golden** (`name2smiles.csv`) — the exact strings the generator emits.
//! * **Property** (everything else) — statements derived from the *name*, not
//!   recorded from the output: carbon counts, which atoms a block may point
//!   at, which constructs are allowed to appear at all. These are the checks
//!   that validate the encoded structure independently.
//!
//! Regenerate the golden file after a deliberate encoding change with
//! `cargo run --example bless_testdata`, then read the diff carefully.

use std::collections::HashSet;
use std::path::PathBuf;

use lipid_notation::{canonicalize, name2smiles, name2structure, smiles2name, smiles_expand};

fn testdata(file: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("testdata")
        .join(file)
}

/// Minimal CSV reader: enough for these files, which hold no quoted fields
/// (SMILES contain commas only inside the `|...|` block, and every such field
/// is quoted — so honour quotes, but nothing fancier).
fn read_csv(file: &str) -> Vec<Vec<String>> {
    let text = std::fs::read_to_string(testdata(file))
        .unwrap_or_else(|e| panic!("read testdata/{file}: {e}"));
    text.lines()
        .skip(1)
        .filter(|l| !l.trim().is_empty())
        .map(|line| {
            let mut fields = Vec::new();
            let mut cur = String::new();
            let mut quoted = false;
            let mut chars = line.chars().peekable();
            while let Some(c) = chars.next() {
                match c {
                    '"' if quoted && chars.peek() == Some(&'"') => {
                        cur.push('"');
                        chars.next();
                    }
                    '"' => quoted = !quoted,
                    ',' if !quoted => fields.push(std::mem::take(&mut cur)),
                    _ => cur.push(c),
                }
            }
            fields.push(cur);
            fields
        })
        .collect()
}

/// Atom symbols in SMILES emission order — an independent reimplementation of
/// the indexing convention every CXSMILES block depends on. Independent on
/// purpose: if the generator's own atom counter drifts, a test sharing that
/// counter would drift with it.
fn atom_symbols(smiles: &str) -> Vec<String> {
    let chars: Vec<char> = smiles.chars().collect();
    let mut out = Vec::new();
    let mut i = 0usize;
    while i < chars.len() {
        match chars[i] {
            '[' => {
                let start = i;
                while i < chars.len() && chars[i] != ']' {
                    i += 1;
                }
                out.push(chars[start..=i.min(chars.len() - 1)].iter().collect());
                i += 1;
            }
            c @ ('B' | 'C' | 'N' | 'O' | 'P' | 'S' | 'F' | 'I' | 'c' | 'n' | 'o' | 's' | 'p'
            | '*') => {
                out.push(c.to_string());
                i += 1;
            }
            _ => i += 1,
        }
    }
    out
}

/// Splits a stored string into base SMILES, the `|...|` blocks, and the
/// trailing size constraint (which lives *outside* the pipes and is not a
/// block).
fn split_stored(stored: &str) -> (&str, &str) {
    match stored.split_once(" |") {
        Some((base, tail)) => (base, tail.split('|').next().unwrap_or("")),
        None => (stored, ""),
    }
}

#[test]
fn golden_strings_still_match() {
    let rows = read_csv("name2smiles.csv");
    assert!(
        rows.len() > 40,
        "corpus looks truncated: {} rows",
        rows.len()
    );
    for row in &rows {
        let (name, cxsmiles, expanded) = (&row[0], &row[1], &row[2]);
        // An empty `cxsmiles` column is a row asserting that the name is
        // *rejected* — a group whose shorthand names something this
        // generator will not invent. Those are worth pinning too: a
        // rejection quietly becoming a guess is the regression that matters
        // most here.
        if cxsmiles.is_empty() {
            assert_eq!(name2smiles(name), None, "{name} should be rejected");
            continue;
        }
        assert_eq!(
            name2smiles(name).as_deref(),
            Some(cxsmiles.as_str()),
            "name2smiles({name})"
        );
        assert_eq!(
            name2structure(name).map(|s| s.smiles).as_deref(),
            Some(expanded.as_str()),
            "name2structure({name})"
        );
        // Expansion and `name2structure` build the same connected molecule.
        assert_eq!(
            smiles_expand(cxsmiles),
            *expanded,
            "{name}: expansion should agree with name2structure"
        );
    }
}

#[test]
fn every_stored_form_canonicalizes_idempotently() {
    for row in read_csv("name2smiles.csv") {
        let (name, stored) = (&row[0], &row[1]);
        if stored.is_empty() {
            continue;
        }
        let canonical = canonicalize(stored)
            .unwrap_or_else(|| panic!("{name}: could not canonicalize {stored}"));
        assert_eq!(
            canonicalize(&canonical).as_deref(),
            Some(canonical.as_str()),
            "{name}: canonical form was not stable"
        );

        // These extensions refer to named Sg variables and atom labels, so
        // atom renumbering must not alter them.
        let trailer = stored.split('|').nth(2).unwrap_or("").trim();
        let canonical_trailer = canonical.split('|').nth(2).unwrap_or("").trim();
        assert_eq!(canonical_trailer, trailer, "{name}: trailer changed");

        let recovered = smiles2name(&canonical)
            .unwrap_or_else(|| panic!("{name}: canonical form should remain readable"));
        let regenerated = name2smiles(&recovered)
            .and_then(|generated| canonicalize(&generated))
            .unwrap_or_else(|| panic!("{name}: recovered name did not regenerate"));
        assert_eq!(
            regenerated, canonical,
            "{name}: recovered name is not canonically equivalent"
        );
    }
}

/// Every chain must contain exactly as many carbons as its name declares, and
/// every mapped index must land on a carbon that no other chain claims.
///
/// The counts in `chains.csv` are read off the names by hand, so this is an
/// independent statement of what the string has to contain — not a recording
/// of what it currently does. It is the check that catches a feature silently
/// dropped from the output.
#[test]
fn chain_carbon_counts_match_the_names() {
    let mut expected: Vec<(String, Vec<(usize, usize)>)> = Vec::new();
    for row in read_csv("chains.csv") {
        let (name, sn, carbons) = (
            row[0].clone(),
            row[1].parse::<usize>().expect("sn"),
            row[2].parse::<usize>().expect("carbons"),
        );
        match expected.last_mut() {
            Some((n, chains)) if *n == name => chains.push((sn, carbons)),
            _ => expected.push((name, vec![(sn, carbons)])),
        }
    }
    assert!(!expected.is_empty(), "chains.csv is empty");

    for (name, chains) in &expected {
        let st = name2structure(name).unwrap_or_else(|| panic!("{name} should resolve"));
        assert!(
            !st.smiles.contains(" |"),
            "{name}: depiction form still carries a CXSMILES block: {}",
            st.smiles
        );
        let symbols = atom_symbols(&st.smiles);
        assert_eq!(symbols.len(), st.atom_count, "{name}: atom_count disagrees");
        assert_eq!(st.chains.len(), chains.len(), "{name}: chain count");

        let mut seen: HashSet<usize> = HashSet::new();
        for (chain, &(sn, carbons)) in st.chains.iter().zip(chains) {
            assert_eq!(chain.sn, sn, "{name}: sn position");
            assert_eq!(chain.carbons.len(), carbons, "{name} sn{sn}: carbon count");
            for (k, &atom) in chain.carbons.iter().enumerate() {
                assert_eq!(
                    symbols.get(atom).map(String::as_str),
                    Some("C"),
                    "{name} sn{sn} C{} maps to atom {atom}, not a carbon, in {}",
                    k + 1,
                    st.smiles
                );
                assert!(
                    seen.insert(atom),
                    "{name}: atom {atom} claimed by two chains"
                );
            }
        }
    }
}

/// A position-variation bond's variable end has to be a `*` dummy carrying
/// exactly one bond, and its candidates have to be chain carbons. Both
/// toolkits enforce the first: CDK ignores an `m:` block that points anywhere
/// else, and RDKit rejects the whole string.
#[test]
fn m_blocks_point_at_a_wildcard_stub() {
    let mut seen = 0;
    for row in read_csv("name2smiles.csv") {
        let (base, blocks) = split_stored(&row[1]);
        let symbols = atom_symbols(base);
        for block in blocks.split(',').filter(|b| b.starts_with("m:")) {
            seen += 1;
            let mut parts = block.split(':').skip(1);
            let idx: usize = parts.next().unwrap().parse().expect("m: atom index");
            assert_eq!(
                symbols.get(idx).map(String::as_str),
                Some("*"),
                "{}: m: targets atom {idx}, not a wildcard, in {base}",
                row[0]
            );
            for site in parts.next().unwrap().split('.') {
                let site: usize = site.parse().expect("m: candidate index");
                assert_eq!(
                    symbols.get(site).map(String::as_str),
                    Some("C"),
                    "{}: m: candidate {site} is not a carbon in {base}",
                    row[0]
                );
            }
        }
    }
    assert!(seen >= 4, "corpus should exercise the m: path");
}

/// Every `snN` label must land on the atom a chain actually hangs from, and
/// the `swappable(...)` token must name exactly the labels present.
///
/// Asserted rather than eyeballed because a misplaced label is invisible: the
/// string still parses, still depicts, and still looks like it says something
/// about sn positions. The labels are also the only in-pipe trace that the
/// trailer has something to say, so an off-by-one here silently unmoors the
/// token from the molecule.
#[test]
fn sn_labels_land_on_the_atom_each_chain_hangs_from() {
    let mut seen = 0;
    for row in read_csv("name2smiles.csv") {
        let (name, stored) = (&row[0], &row[1]);
        if stored.is_empty() {
            continue;
        }
        let (base, blocks) = split_stored(stored);
        let trailer = stored.split('|').next_back().unwrap_or("").trim();
        let swappable = trailer.contains("swappable(");

        if !swappable {
            assert!(
                !blocks.contains('$'),
                "{name}: sn labels with no swappable token"
            );
            continue;
        }
        seen += 1;

        let labels = blocks
            .split_once('$')
            .and_then(|(_, r)| r.split_once('$'))
            .map(|(l, _)| l)
            .unwrap_or_else(|| panic!("{name}: swappable token with no $...$ labels"));

        let symbols = atom_symbols(base);
        let mut found: Vec<&str> = Vec::new();
        for (i, slot) in labels.split(';').enumerate() {
            if slot.is_empty() {
                continue;
            }
            assert!(
                slot.starts_with("sn") && slot[2..].parse::<u32>().is_ok(),
                "{name}: unexpected atom label {slot}"
            );
            // A chain hangs from its ester/ether oxygen on every class that
            // has an sn position to be ambiguous about.
            assert_eq!(
                symbols.get(i).map(String::as_str),
                Some("O"),
                "{name}: {slot} labels atom {i}, which is not a linking oxygen, in {base}"
            );
            found.push(slot);
        }

        let named: Vec<&str> = trailer
            .split("swappable(")
            .nth(1)
            .and_then(|t| t.split(')').next())
            .unwrap_or("")
            .split(',')
            .collect();
        found.sort();
        let mut named_sorted = named.clone();
        named_sorted.sort();
        assert_eq!(
            found, named_sorted,
            "{name}: swappable names {named:?} but the string labels {found:?}"
        );
        assert!(found.len() >= 2, "{name}: nothing to swap");
    }
    assert!(seen >= 4, "corpus should exercise the swappable path");
}

/// A fully determined name must come out as plain SMILES with no tail at all,
/// because the presence of a tail is the signal that something was
/// undetermined.
#[test]
fn determined_names_carry_no_cxsmiles_tail() {
    for name in [
        "FA 18:0",
        "FA 18:1(9Z)",
        "FA 18:1(9)",
        "FA 18:0;5OH",
        "PC 16:0/18:1(9)",
        "TG 16:0/18:1(9)/18:2(9,12)",
    ] {
        let s = name2smiles(name).unwrap_or_else(|| panic!("{name} should resolve"));
        assert!(!s.contains(" |"), "{name} should be plain SMILES: {s}");
    }
}

/// Multi-chain shorthand has many chain realizations; picking one would be a
/// fabrication, so nothing is emitted.
#[test]
fn multi_chain_shorthand_is_rejected() {
    for name in ["PC 34:1", "PC 19:2", "DG 36:2", "TG 54:3"] {
        assert_eq!(name2smiles(name), None, "{name} should be rejected");
    }
    // Single-chain classes have nothing to be ambiguous between.
    assert!(name2smiles("FA 18:1").is_some());
}

/// The illustrative strings used by the documentation.
#[test]
fn doc_examples_still_match() {
    for row in read_csv("doc_examples.csv") {
        let (name, cxsmiles) = (&row[0], &row[1]);
        assert_eq!(
            name2smiles(name).as_deref(),
            Some(cxsmiles.as_str()),
            "doc example {name}"
        );
    }
}

/// The circular test: every name in the corpus, through SMILES and back.
///
/// This is the check `smiles2name` exists to make possible, and it is a
/// different kind of check from everything above it. The golden strings pin
/// what the generator emits; this pins that what it emits still *means* the
/// name it came from. A miscounted chain, a double bond written one carbon
/// off, an `Sg:` run whose constraint no longer matches its markers — none of
/// those change whether a string parses, and a literal comparison against a
/// recorded expectation cannot see any of them. A round trip can.
///
/// Two properties are asserted, and the weaker one is the important one:
///
/// 1. **Every name comes back** — no corpus name reads as unrecognizable.
/// 2. **The name that comes back regenerates the same string.** This holds
///    even where the first property's exact-equality version does not, and it
///    is what makes the round trip a validation rather than a formatting
///    check.
#[test]
fn every_name_survives_the_round_trip() {
    let mut inexact = Vec::new();

    for row in read_csv("name2smiles.csv") {
        let (name, cxsmiles) = (&row[0], &row[1]);
        if cxsmiles.is_empty() {
            continue; // a row pinning a rejection; there is nothing to invert
        }

        let recovered = smiles2name(cxsmiles)
            .unwrap_or_else(|| panic!("{name}: no name recovered from {cxsmiles}"));

        // Compared canonically, not textually: the reverse direction works in
        // canonical space, so it may return a different *spelling* of the same
        // structure. `FA 18:1(6Z);[14-18cy5:1(15)]` comes back as `…(17)`
        // because a cyclopentene ring is symmetric about the carbon it hangs
        // from — same molecule, other name for it.
        assert_eq!(
            name2smiles(&recovered).and_then(|s| canonicalize(&s)),
            canonicalize(cxsmiles),
            "{name}: recovered as {recovered}, which is a different structure"
        );

        // A second trip must change nothing: whatever the first one
        // canonicalized is already canonical.
        let again = smiles2name(cxsmiles).expect("still recoverable");
        assert_eq!(recovered, again, "{name}: round trip is not a fixed point");

        if recovered != *name {
            inexact.push(format!("{name} -> {recovered}"));
        }
    }

    // The names that do *not* come back verbatim are the places where
    // `name2smiles` is genuinely not injective. Each is a documented
    // limitation rather than a bug, so the set is pinned: a new entry here
    // means information started being lost somewhere it previously was not.
    inexact.sort();
    assert_eq!(
        inexact,
        [
            // A cyclopentene ring is symmetric about the carbon it hangs
            // from, so a double bond at 15-16 and one at 17-18 are the same
            // molecule under two names.
            "FA 18:1(6Z);[14-18cy5:1(15)] -> FA 18:1(6Z);[14-18cy5:1(17)]",
            // `cy5` and `cy5:0` are the same ring; the count is written out.
            "FA 20:2(5,8);[11-15cy5;13OH];18OH -> FA 20:2(5,8);[11-15cy5:0;13OH];18OH",
            // sn-1 and sn-3 are a *stereochemical* distinction, and this
            // generator writes no stereocenter at glycerol C2 — so the two
            // primary positions are indistinguishable in the structure and
            // the reverse direction cannot tell which chain was which. Both
            // names produce the same molecule. Worth knowing: a `/`-joined
            // TG or DG asserts an sn-1/sn-3 order its own SMILES does not
            // actually encode.
            "TG 16:0/18:0;5Me/18:1(9) -> TG 18:1(9)/18:0;5Me/16:0",
        ],
        "the set of names that lose information changed"
    );

    // Unresolved `_` assignments survive verbatim, including alongside Sg/m.
    for name in [
        "PC 16:0_18:1",
        "PC 16:0_18:1(9);OH",
        "DG 16:0_18:1",
        "DG 18:1(9);5OH_18:1;OH",
        "TG 16:0_18:1_18:2",
        "TG 18:0_18:1_18:2",
    ] {
        let smi = name2smiles(name).unwrap_or_else(|| panic!("{name} should resolve"));
        assert_eq!(smiles2name(&smi).as_deref(), Some(name), "{name}");
    }
}

/// The expanded depiction strings are ordinary SMILES with every `Sg:` run
/// resolved to one concrete split, so they must read back as a *fully
/// localized* name — one that regenerates the expansion rather than the
/// ambiguous original. Getting the original back would mean the expansion had
/// smuggled its arbitrary choice in as a determination.
#[test]
fn expanded_depictions_read_back_as_localized_names() {
    let mut checked = 0;
    for row in read_csv("name2smiles.csv") {
        let (cxsmiles, expanded) = (&row[1], &row[2]);
        if cxsmiles.is_empty() || !cxsmiles.contains("Sg:") {
            continue;
        }
        let Some(recovered) = smiles2name(expanded) else {
            continue; // expansion can leave shapes the templates do not cover
        };
        checked += 1;
        assert_eq!(
            name2smiles(&recovered).as_deref(),
            Some(expanded.as_str()),
            "{}: {recovered} does not regenerate its own expansion",
            row[0]
        );
        assert!(
            !recovered.contains("Sg:"),
            "{}: expansion should not read back as ambiguous",
            row[0]
        );
    }
    assert!(checked >= 5, "corpus should exercise the expanded path");
}

/// Every generator string quoted in the README must match current output.
///
/// A tail containing `…` is a deliberate elision and is skipped.
#[test]
fn doc_prose_quotes_current_output() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let by_base: Vec<(String, String, String)> = read_csv("name2smiles.csv")
        .iter()
        .filter(|row| !row[1].is_empty())
        .map(|row| {
            let base = row[1].split(" |").next().unwrap_or("").to_string();
            (base, row[0].clone(), row[1].clone())
        })
        .filter(|(base, _, _)| base.len() >= 25)
        .collect();

    let mut stale = Vec::new();
    for doc in ["README.md"] {
        let Ok(text) = std::fs::read_to_string(root.join(doc)) else {
            continue;
        };
        for (lineno, line) in text.lines().enumerate() {
            for (base, name, _) in &by_base {
                let Some(after) = line.split(base.as_str()).nth(1) else {
                    continue;
                };
                // Only compare where the doc actually shows the tail, and
                // never across an elision.
                if !after.starts_with(" |") || after.contains('…') {
                    continue;
                }
                // Trim the punctuation that wraps a quotation, but not the
                // `)` that closes a `constrain(...)` token.
                // An inline annotation is set off by a run of spaces.
                let mut quoted = after.split("  ").next().unwrap_or(after).trim();
                loop {
                    let trimmed = match quoted.chars().next_back() {
                        Some('"' | '`') => &quoted[..quoted.len() - 1],
                        Some(')') if quoted.matches(')').count() > quoted.matches('(').count() => {
                            &quoted[..quoted.len() - 1]
                        }
                        _ => break,
                    };
                    quoted = trimmed.trim_end();
                }
                // Several names share a base SMILES and differ only in their
                // blocks — `PC 16:0/18:1` and `PC 16:0_18:1` now do. Any of
                // them matching is a match.
                let mut expected: Vec<&str> = by_base
                    .iter()
                    .filter(|(b, _, _)| b == base)
                    .map(|(b, _, cur)| cur.strip_prefix(b.as_str()).unwrap_or("").trim())
                    .collect();
                expected.sort();
                expected.dedup();
                if !expected.contains(&quoted) {
                    stale.push(format!(
                        "{doc}:{}: {name}\n     doc: {quoted}\n    real: {}",
                        lineno + 1,
                        expected.join("\n        or: ")
                    ));
                }
            }
        }
    }
    assert!(
        stale.is_empty(),
        "documentation quotes output the generator no longer produces:\n{}",
        stale.join("\n")
    );
}

// ---------- corpus regeneration ----------

/// Quote a CSV field iff it needs it.
fn csv_field(s: &str) -> String {
    if s.contains(',') || s.contains('"') || s.contains('\n') {
        format!("\"{}\"", s.replace('"', "\"\""))
    } else {
        s.to_string()
    }
}

/// Notes worth preserving across a rewrite — they say *why* a row is in the
/// file, which the generator cannot regenerate.
fn corpus_note(name: &str) -> &'static str {
    match name {
        "FA 18:0;ep(5)" => "accepted ;ep(pos) alias for the Table 1A Ep group",
        _ => "",
    }
}

/// Rewrites the generated columns of `name2smiles.csv` and `doc_examples.csv`
/// from the names already listed there. Run after a *deliberate* encoding
/// change, then read the diff — that diff is the only review a golden file
/// ever gets:
///
/// ```text
/// cargo test --test testdata -- --ignored
/// git diff testdata/
/// ```
///
/// The `name` column is the corpus definition and is never rewritten; add a
/// case by adding a row with the generated columns left empty. The
/// hand-written counts in `chains.csv` are never touched — they exist
/// precisely to be independent of what the code emits.
#[test]
#[ignore = "rewrites the checked-in testdata corpus"]
fn write_testdata_corpus() {
    for (file, with_note) in [("name2smiles.csv", false), ("doc_examples.csv", true)] {
        let path = testdata(file);
        let mut out = String::from(if with_note {
            "name,cxsmiles,expanded,note\n"
        } else {
            "name,cxsmiles,expanded\n"
        });

        for row in read_csv(file) {
            let name = &row[0];
            let cxsmiles = name2smiles(name).unwrap_or_default();
            let expanded = name2structure(name).map(|s| s.smiles).unwrap_or_default();
            if cxsmiles.is_empty() {
                // Not fatal: a row may exist to pin that a name is rejected.
                println!("  note: {name} produced no structure");
            }
            out.push_str(&csv_field(name));
            out.push(',');
            out.push_str(&csv_field(&cxsmiles));
            out.push(',');
            out.push_str(&csv_field(&expanded));
            if with_note {
                out.push(',');
                out.push_str(&csv_field(corpus_note(name)));
            }
            out.push('\n');
        }

        std::fs::write(&path, out).unwrap_or_else(|e| panic!("write {path:?}: {e}"));
        println!("rewrote testdata/{file}");
    }
    println!("\nnow read the diff: git diff testdata/");
}
