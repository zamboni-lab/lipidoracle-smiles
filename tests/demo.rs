//! Generates the checked-in `demo.html` gallery.
//!
//! Run `cargo test --test demo -- --ignored` after changing examples or
//! conversion behaviour. The ordinary test validates the examples without
//! modifying the working tree.

use std::path::PathBuf;

use lipid_notation::{canonicalize, name2smiles, smiles2name, smiles_for_depiction};

const CDK_DEPICT_SVG: &str = "https://www.simolecule.com/cdkdepict/depict/bow/svg?smi=";

struct Example {
    group: &'static str,
    shorthand: &'static str,
    valid: bool,
}

const EXAMPLES: &[Example] = &[
    Example {
        group: "Fully determined structures",
        shorthand: "FA 18:1(9Z)",
        valid: true,
    },
    Example {
        group: "Fully determined structures",
        shorthand: "FA 20:4(5Z,8Z,12E,14Z);11OH",
        valid: true,
    },
    Example {
        group: "Fully determined structures",
        shorthand: "PC 16:0/18:1(9Z)",
        valid: true,
    },
    Example {
        group: "Fully determined structures",
        shorthand: "TG 16:0/18:0;5Me/18:1(9)",
        valid: true,
    },
    Example {
        group: "Fully determined structures",
        shorthand: "TG 16:0/18:1(9)/18:2(9,12)",
        valid: true,
    },
    Example {
        group: "Unlocalized lipid features",
        shorthand: "FA 20:4",
        valid: true,
    },
    Example {
        group: "Unlocalized lipid features",
        shorthand: "FA 20:4(5,8,12,14);OH",
        valid: true,
    },
    Example {
        group: "Unlocalized lipid features",
        shorthand: "PC 16:0_18:1",
        valid: true,
    },
    Example {
        group: "Unlocalized lipid features",
        shorthand: "PC O-16:1_18:2;OH",
        valid: true,
    },
    Example {
        group: "Combined rings, groups, and uncertainty",
        shorthand: "FA 20:3(5,8,11);(OH)2",
        valid: true,
    },
    Example {
        group: "Combined rings, groups, and uncertainty",
        shorthand: "FA 20:2(5,8);[11-15cy5;13OH];OH",
        valid: true,
    },
    Example {
        group: "Combined rings, groups, and uncertainty",
        shorthand: "PC 16:1_18:2;9Ep;OH",
        valid: true,
    },
    Example {
        group: "Combined rings, groups, and uncertainty",
        shorthand: "TG 18:1(9);5OH_18:2;9Ep_18:1",
        valid: true,
    },
    Example {
        group: "Combined rings, groups, and uncertainty",
        shorthand: "Cer d18:1(4)/16:1;OH",
        valid: true,
    },
    Example {
        group: "Functional-group and ring notation",
        shorthand: "FA 16:0;3Me,7Me,11Me,15Me",
        valid: true,
    },
    Example {
        group: "Functional-group and ring notation",
        shorthand: "FA 18:1(9);12NO2",
        valid: true,
    },
    Example {
        group: "Functional-group and ring notation",
        shorthand: "FA 18:0;9Ep",
        valid: true,
    },
    Example {
        group: "Functional-group and ring notation",
        shorthand: "FA 19:0;[11-13cy3:0]",
        valid: true,
    },
    Example {
        group: "Functional-group and ring notation",
        shorthand: "FA 19:0;[9-11cy3:1(9)]",
        valid: true,
    },
    Example {
        group: "Functional-group and ring notation",
        shorthand: "FA 20:2(5,8);[11-15cy5;13OH];18OH",
        valid: true,
    },
    Example {
        group: "Headgroup diversity",
        shorthand: "PC P-16:0/18:1(9)",
        valid: true,
    },
    Example {
        group: "Headgroup diversity",
        shorthand: "CAR 18:1(9)",
        valid: true,
    },
    Example {
        group: "Headgroup diversity",
        shorthand: "NAE 20:4(5,8,11,14)",
        valid: true,
    },
    Example {
        group: "Headgroup diversity",
        shorthand: "AMP-FA 20:4(5,8,11,14);15OH",
        valid: true,
    },
    Example {
        group: "Headgroup diversity",
        shorthand: "Cer d18:1(4)/16:0",
        valid: true,
    },
    Example {
        group: "Headgroup diversity",
        shorthand: "SM d18:1(4)/16:0",
        valid: true,
    },
    Example {
        group: "Headgroup diversity",
        shorthand: "HexCer d18:1(4)/16:0",
        valid: true,
    },
    Example {
        group: "Headgroup diversity",
        shorthand: "S1P d18:1(4)",
        valid: true,
    },
    Example {
        group: "Intentionally unsupported shorthand",
        shorthand: "PC 34:1",
        valid: false,
    },
    Example {
        group: "Intentionally unsupported shorthand",
        shorthand: "FA 18:1(9);O2",
        valid: false,
    },
    Example {
        group: "Intentionally unsupported shorthand",
        shorthand: "FA 18:0;1OMe",
        valid: false,
    },
    Example {
        group: "Intentionally unsupported shorthand",
        shorthand: "Gal-Glc-Cer d18:1(4)/16:0",
        valid: false,
    },
];

#[derive(Debug)]
struct ResultRow {
    generated: Option<String>,
    canonical: Option<String>,
    recovered: Option<String>,
    direct_depiction: Option<String>,
    prepared_depiction: Option<String>,
}

fn evaluate(example: &Example) -> ResultRow {
    let generated = name2smiles(example.shorthand);
    let canonical = generated.as_deref().and_then(canonicalize);
    let recovered = canonical.as_deref().and_then(smiles2name);
    let direct_depiction = generated.clone();
    let prepared_depiction = generated.as_deref().map(smiles_for_depiction);

    ResultRow {
        generated,
        canonical,
        recovered,
        direct_depiction,
        prepared_depiction,
    }
}

fn validate(example: &Example, row: &ResultRow) {
    if !example.valid {
        assert!(
            row.generated.is_none(),
            "{} should be rejected",
            example.shorthand
        );
        assert!(row.canonical.is_none());
        assert!(row.recovered.is_none());
        assert!(row.direct_depiction.is_none());
        assert!(row.prepared_depiction.is_none());
        return;
    }

    let generated = row
        .generated
        .as_deref()
        .unwrap_or_else(|| panic!("{} should generate", example.shorthand));
    let canonical = row
        .canonical
        .as_deref()
        .unwrap_or_else(|| panic!("{} should canonicalize", example.shorthand));
    let recovered = row
        .recovered
        .as_deref()
        .unwrap_or_else(|| panic!("{} should read back", example.shorthand));
    let regenerated = name2smiles(recovered)
        .as_deref()
        .and_then(canonicalize)
        .expect("recovered name should regenerate");

    assert_eq!(
        regenerated, canonical,
        "{} recovered as {recovered}, which changed the structure",
        example.shorthand
    );
    assert!(row.direct_depiction.is_some());
    assert!(row.prepared_depiction.is_some());
    assert!(generated.starts_with(|c: char| c.is_ascii_alphabetic() || c == '[' || c == '*'));
}

/// Whether the name that came back is the name that went in.
///
/// A mismatch is not a failure. `name2smiles` is not injective, so more than
/// one name can describe the same structure, and the reverse direction returns
/// whichever spelling it derives from the molecule. What matters is that the
/// recovered name regenerates the same structure — the round-trip test in
/// `tests/testdata.rs` asserts exactly that, and this column shows where the
/// spelling nevertheless differs.
fn name_identity(shorthand: &str, recovered: Option<&str>) -> String {
    match recovered {
        Some(name) if name == shorthand => "yes".to_string(),
        Some(name) => format!("as {}", html_code(Some(name))),
        None => "&mdash;".to_string(),
    }
}

fn html_code(value: Option<&str>) -> String {
    let value = value.unwrap_or("—");
    format!(
        "<code>{}</code>",
        value
            .replace('&', "&amp;")
            .replace('<', "&lt;")
            .replace('>', "&gt;")
            .replace('|', "&#124;")
    )
}

fn percent_encode(value: &str) -> String {
    let mut out = String::with_capacity(value.len() * 3);
    for byte in value.bytes() {
        if byte.is_ascii_alphanumeric() || matches!(byte, b'-' | b'_' | b'.' | b'~') {
            out.push(byte as char);
        } else {
            out.push('%');
            out.push_str(&format!("{byte:02X}"));
        }
    }
    out
}

fn cdk_depiction(smiles: Option<&str>) -> String {
    let Some(smiles) = smiles else {
        return "—".to_string();
    };
    let url = format!(
        "{}{encoded}&amp;showtitle=false",
        CDK_DEPICT_SVG,
        encoded = percent_encode(smiles)
    );
    format!("<img src=\"{url}\" alt=\"CDK depiction\" />")
}

fn render_demo() -> String {
    let mut out = String::from(
        r#"<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>lipid_notation demo</title>
<style>
  :root { color-scheme: light dark; font-family: system-ui, sans-serif; }
  body { margin: 2rem; }
  .table-wrap { overflow-x: auto; border: 1px solid #8885; }
  table { border-collapse: collapse; min-width: 2560px; width: 100%; }
  th, td { border: 1px solid #8885; padding: .7rem; vertical-align: top; text-align: left; }
  th { position: sticky; top: 0; background: Canvas; }
  code { white-space: pre-wrap; overflow-wrap: anywhere; font-size: .8rem; }
  .depiction { min-width: 660px; text-align: center; }
  .depiction img { width: 640px; max-width: none; height: auto; background: white; }
  .identity { text-align: center; white-space: nowrap; }
</style>
</head>
<body>
<h1>lipid_notation demo</h1>
<p>This gallery is generated by <code>tests/demo.rs</code>. Regenerate it with
<code>cargo test --test demo -- --ignored</code>.</p>
<p>Each valid shorthand is converted with <code>name2smiles</code>, canonicalized
with <code>canonicalize</code>, and read back with <code>smiles2name</code>.
The two SVG columns compare direct rendering of the CXSMILES with the output of
<code>smiles_for_depiction</code>. The latter canonicalizes and reindexes the CXSMILES,
then places each <code>m:</code> group over the nearest unused side-chain single bond.
The full uncertainty encoding remains in the data columns.</p>
"#,
    );

    let mut group = "";
    for example in EXAMPLES {
        if example.group != group {
            if !group.is_empty() {
                out.push_str("</tbody></table></div></section>\n");
            }
            group = example.group;
            out.push_str(&format!(
                "<section><h2>{group}</h2><div class=\"table-wrap\"><table>\n"
            ));
            out.push_str("<thead><tr><th>Test shorthand</th><th><code>name2smiles</code> output</th><th>Canonicalized CXSMILES</th><th><code>smiles2name</code> result</th><th>Name recovered?</th><th>Direct CXSMILES depiction</th><th><code>smiles_for_depiction</code> depiction</th></tr></thead><tbody>\n");
        }

        let row = evaluate(example);
        validate(example, &row);
        out.push_str(&format!(
            "<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td class=\"identity\">{}</td><td class=\"depiction\">{}</td><td class=\"depiction\">{}</td></tr>\n",
            html_code(Some(example.shorthand)),
            html_code(row.generated.as_deref()),
            html_code(row.canonical.as_deref()),
            html_code(row.recovered.as_deref()),
            name_identity(example.shorthand, row.recovered.as_deref()),
            cdk_depiction(row.direct_depiction.as_deref()),
            cdk_depiction(row.prepared_depiction.as_deref()),
        ));
    }
    out.push_str("</tbody></table></div></section>\n</body>\n</html>\n");
    out
}

#[test]
fn representative_examples_are_convertible_or_explicitly_rejected() {
    for example in EXAMPLES {
        validate(example, &evaluate(example));
    }
}

#[test]
fn checked_in_demo_matches_the_generator() {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("demo.html");
    let checked_in = std::fs::read_to_string(&path)
        .unwrap_or_else(|error| panic!("read {}: {error}", path.display()));
    assert_eq!(
        checked_in,
        render_demo(),
        "run `cargo test --test demo -- --ignored`"
    );
}

#[test]
#[ignore = "writes the checked-in demo.html gallery"]
fn write_demo_html() {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("demo.html");
    std::fs::write(&path, render_demo()).unwrap_or_else(|error| {
        panic!("write {}: {error}", path.display());
    });
}
