//! Regenerates the golden columns of `testdata/name2smiles.csv` and
//! `testdata/doc_examples.csv` from the names already listed there.
//!
//! Run after a *deliberate* encoding change, then read the diff carefully —
//! that diff is the only review a golden file ever gets:
//!
//! ```text
//! cargo run --example bless_testdata
//! git diff testdata/
//! ```
//!
//! The `name` column is the corpus definition and is never rewritten; add new
//! cases by adding a row with an empty `cxsmiles`/`expanded` and blessing. The
//! hand-written expectations in `testdata/chains.csv` are never touched by
//! this: they exist precisely to be independent of what the code emits.

use std::path::{Path, PathBuf};

use lipid_notation::{name2smiles, name2structure};

fn testdata(file: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("testdata")
        .join(file)
}

/// Quote a CSV field iff it needs it.
fn field(s: &str) -> String {
    if s.contains(',') || s.contains('"') || s.contains('\n') {
        format!("\"{}\"", s.replace('"', "\"\""))
    } else {
        s.to_string()
    }
}

/// The first column of every data row, which is the corpus definition.
fn names(path: &Path) -> Vec<String> {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("read {path:?}: {e}"));
    text.lines()
        .skip(1)
        .filter(|l| !l.trim().is_empty())
        .map(|line| {
            if let Some(rest) = line.strip_prefix('"') {
                rest.split_once('"')
                    .map(|(n, _)| n)
                    .unwrap_or(rest)
                    .to_string()
            } else {
                line.split(',').next().unwrap_or("").to_string()
            }
        })
        .collect()
}

/// Notes worth preserving across a bless — they describe *why* a row is in the
/// file, which the generator cannot regenerate.
fn note(name: &str) -> &'static str {
    match name {
        "FA 18:0;ep(5)" => "accepted ;ep(pos) alias for the Table 1A Ep group",
        _ => "",
    }
}

fn main() {
    for (file, with_note) in [("name2smiles.csv", false), ("doc_examples.csv", true)] {
        let path = testdata(file);
        let mut out = if with_note {
            String::from("name,cxsmiles,expanded,note\n")
        } else {
            String::from("name,cxsmiles,expanded\n")
        };
        let mut missing = Vec::new();

        for name in names(&path) {
            let cxsmiles = name2smiles(&name).unwrap_or_default();
            let expanded = name2structure(&name).map(|s| s.smiles).unwrap_or_default();
            if cxsmiles.is_empty() {
                missing.push(name.clone());
            }
            out.push_str(&field(&name));
            out.push(',');
            out.push_str(&field(&cxsmiles));
            out.push(',');
            out.push_str(&field(&expanded));
            if with_note {
                out.push(',');
                out.push_str(&field(note(&name)));
            }
            out.push('\n');
        }

        std::fs::write(&path, out).unwrap_or_else(|e| panic!("write {path:?}: {e}"));
        println!("blessed testdata/{file}");
        for name in missing {
            // Not fatal: a row may exist to document that a name is rejected.
            println!("  note: {name} produced no structure");
        }
    }
    println!("\nnow read the diff: git diff testdata/");
}
