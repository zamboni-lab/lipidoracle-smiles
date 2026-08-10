//! Corpus tests driven by the CSVs in `testdata/`.
//!
//! Two kinds of check live here, and the distinction matters:
//!
//! * **Golden** (`name2smiles.csv`) — the exact strings the generator emits.
//!   Catches regressions. Cannot catch a misconception: if an expected string
//!   was wrong from the start, comparing against it passes forever.
//! * **Property** (everything else) — statements derived from the *name*, not
//!   recorded from the output: carbon counts, which atoms a block may point
//!   at, which constructs are allowed to appear at all. These are the checks
//!   that catch a block that was wrong from day one, and three of them were.
//!
//! Regenerate the golden file after a deliberate encoding change with
//! `cargo run --example bless_testdata`, then read the diff carefully.

use std::collections::HashSet;
use std::path::PathBuf;

use shorthand2smiles::{expand_cxsmiles_for_depiction, name2smiles, name2structure};

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
            c @ ('C' | 'N' | 'O' | 'P' | 'S' | 'F' | 'I' | 'c' | 'n' | 'o' | 's' | 'p' | '*') => {
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
        // `expand_cxsmiles_for_depiction` resolves `Sg:` runs but cannot
        // resolve an `RG:` scheme — it has no basis for choosing which
        // alternative goes where, so it strips the block and leaves the bare
        // backbone. `name2structure` builds the resolved layout instead. The
        // two therefore agree exactly when no `RG:` block is present, and
        // that boundary is worth pinning down.
        let depicted = expand_cxsmiles_for_depiction(cxsmiles);
        if cxsmiles.contains("RG:") {
            assert!(
                depicted.contains('*') && !depicted.contains(" |"),
                "{name}: expansion should strip the RG: block and leave the slots: {depicted}"
            );
        } else {
            assert_eq!(
                depicted, *expanded,
                "{name}: expansion should agree with name2structure when there is no RG: block"
            );
        }
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

/// Every `R1` slot of the `$...$` label block must land on a `*` atom, and
/// there must be one `RG:` definition per interchangeable chain.
///
/// Asserted rather than eyeballed because CDK renders a misplaced label
/// without complaint — producing a picture that looks right and means
/// something else.
#[test]
fn r_group_labels_land_on_wildcard_atoms() {
    let mut seen_rg = 0;
    for row in read_csv("name2smiles.csv") {
        let (name, stored) = (&row[0], &row[1]);
        let (base, blocks) = split_stored(stored);
        if !blocks.contains("RG:") {
            assert!(
                !blocks.contains('$'),
                "{name}: label block with no RG: definitions"
            );
            continue;
        }
        seen_rg += 1;
        let symbols = atom_symbols(base);
        let labels = blocks
            .split_once('$')
            .and_then(|(_, r)| r.split_once('$'))
            .map(|(l, _)| l)
            .unwrap_or_else(|| panic!("{name}: RG: block with no $...$ labels"));

        let mut sites = 0;
        for (i, slot) in labels.split(';').enumerate() {
            if slot.is_empty() {
                continue;
            }
            assert_eq!(slot, "R1", "{name}: unexpected atom label {slot}");
            assert_eq!(
                symbols.get(i).map(String::as_str),
                Some("*"),
                "{name}: R1 label at slot {i} lands on {:?}, not a wildcard, in {base}",
                symbols.get(i)
            );
            sites += 1;
        }
        let defs = blocks.matches("},{").count() + 1;
        assert_eq!(sites, defs, "{name}: one R1 site per RG: definition");
        assert_eq!(
            sites,
            base.matches('*').count(),
            "{name}: every wildcard in the backbone should be labelled"
        );
    }
    assert!(seen_rg >= 4, "corpus should exercise the RG: path");
}

/// `ctu:` is a ChemAxon *query* feature for matching any configuration, and
/// `f:` groups components into one entity. Neither belongs in a structure:
/// an undecorated `C=C` is already unspecified geometry, and `f:` expresses
/// *and*, never the *or* that unresolved sn-regiochemistry needs. Both were
/// emitted by earlier revisions; see `dev/SMILES.md` §3.
#[test]
fn no_block_is_a_query_or_a_fragment_group() {
    for row in read_csv("name2smiles.csv") {
        let (name, stored) = (&row[0], &row[1]);
        let (_, blocks) = split_stored(stored);
        assert!(
            !blocks.contains("ctu:"),
            "{name} still emits ctu:: {stored}"
        );
        assert!(!blocks.contains("f:"), "{name} still emits f:: {stored}");
    }
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

/// The illustrative strings quoted in `dev/`, including one that documents a
/// known bug rather than intended behaviour.
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
