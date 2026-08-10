//! Generates the checked-in `demo.html` gallery.
//!
//! Run `cargo test --test demo -- --ignored` after changing examples or
//! conversion behaviour. The ordinary test validates the examples without
//! modifying the working tree.

use std::path::PathBuf;

use lipid_notation::{canonicalize_cxsmiles, name2smiles, smiles2name, smiles_for_depiction};

const CDK_DEPICT_SVG: &str = "https://www.simolecule.com/cdkdepict/depict/bow/svg?smi=";

struct Example {
    group: &'static str,
    shorthand: &'static str,
    note: &'static str,
    valid: bool,
}

const EXAMPLES: &[Example] = &[
    Example {
        group: "Fully determined structures",
        shorthand: "FA 18:1(9Z)",
        note: "Oleate-style monounsaturated fatty acid; Z geometry is explicit.",
        valid: true,
    },
    Example {
        group: "Fully determined structures",
        shorthand: "FA 20:4(5Z,8Z,12E,14Z);11OH",
        note: "Four localized double bonds with mixed geometry and one localized hydroxyl.",
        valid: true,
    },
    Example {
        group: "Fully determined structures",
        shorthand: "PC 16:0/18:1(9Z)",
        note: "Resolved sn assignment (`/`) and double-bond geometry.",
        valid: true,
    },
    Example {
        group: "Fully determined structures",
        shorthand: "TG 16:0/18:1(9)/18:2(9,12)",
        note: "Triacylglycerol with three explicit, ordered chains.",
        valid: true,
    },
    Example {
        group: "Unlocalized lipid features",
        shorthand: "FA 20:4",
        note: "Four double bonds are declared but not localized; `Sg:` and `constrain(...)` retain that uncertainty.",
        valid: true,
    },
    Example {
        group: "Unlocalized lipid features",
        shorthand: "FA 20:4(5,8,12,14);OH",
        note: "Localized double bonds with a hydroxyl whose position is represented by `m:`.",
        valid: true,
    },
    Example {
        group: "Unlocalized lipid features",
        shorthand: "PC 16:0_18:1",
        note: "Known chains with unresolved sn assignment (`_`); atom labels and `swappable(...)` preserve it.",
        valid: true,
    },
    Example {
        group: "Unlocalized lipid features",
        shorthand: "PC O-16:1_18:2;OH",
        note: "Ether-linked phosphatidylcholine combining `Sg:`, `m:`, and unresolved sn assignment.",
        valid: true,
    },
    Example {
        group: "Combined rings, groups, and uncertainty",
        shorthand: "FA 20:3(5,8,11);(OH)2",
        note: "Three localized double bonds plus two independently unlocalized hydroxyl groups (`m:` twice).",
        valid: true,
    },
    Example {
        group: "Combined rings, groups, and uncertainty",
        shorthand: "FA 20:2(5,8);[11-15cy5;13OH];OH",
        note: "A cyclopentane and localized in-ring hydroxyl combined with an unlocalized hydroxyl (`m:`).",
        valid: true,
    },
    Example {
        group: "Combined rings, groups, and uncertainty",
        shorthand: "PC 16:1_18:2;9Ep;OH",
        note: "Unresolved sn assignment with `Sg:` regions on both chains, an epoxide, and a position-variable hydroxyl.",
        valid: true,
    },
    Example {
        group: "Combined rings, groups, and uncertainty",
        shorthand: "TG 18:1(9);5OH_18:2;9Ep_18:1",
        note: "Three unresolved glycerolipid chains combining a localized hydroxyl, epoxide, and unlocalized double-bond regions.",
        valid: true,
    },
    Example {
        group: "Combined rings, groups, and uncertainty",
        shorthand: "Cer d18:1(4)/16:1;OH",
        note: "Ceramide with a localized sphingoid double bond and an N-acyl chain carrying both `Sg:` and `m:` ambiguity.",
        valid: true,
    },
    Example {
        group: "Functional-group and ring notation",
        shorthand: "FA 16:0;3Me,7Me,11Me,15Me",
        note: "Phytanic-acid style methyl branching using Table 1A substituents.",
        valid: true,
    },
    Example {
        group: "Functional-group and ring notation",
        shorthand: "FA 18:1(9);12NO2",
        note: "Nitro substituent on a localized unsaturated chain.",
        valid: true,
    },
    Example {
        group: "Functional-group and ring notation",
        shorthand: "FA 18:0;9Ep",
        note: "Epoxide ring across adjacent chain carbons.",
        valid: true,
    },
    Example {
        group: "Functional-group and ring notation",
        shorthand: "FA 19:0;[11-13cy3:0]",
        note: "Cyclopropane ring notation (lactobacillic-acid style).",
        valid: true,
    },
    Example {
        group: "Functional-group and ring notation",
        shorthand: "FA 19:0;[9-11cy3:1(9)]",
        note: "Cyclopropene ring with its own localized double bond (sterculic-acid style).",
        valid: true,
    },
    Example {
        group: "Functional-group and ring notation",
        shorthand: "FA 20:2(5,8);[11-15cy5;13OH];18OH",
        note: "Cyclopentane with an in-ring hydroxyl plus a second chain hydroxyl.",
        valid: true,
    },
    Example {
        group: "Headgroup diversity",
        shorthand: "PC P-16:0/18:1(9)",
        note: "Plasmalogen phosphatidylcholine with a vinyl-ether chain.",
        valid: true,
    },
    Example {
        group: "Headgroup diversity",
        shorthand: "CAR 18:1(9)",
        note: "Acylcarnitine ester with a zwitterionic carnitine headgroup.",
        valid: true,
    },
    Example {
        group: "Headgroup diversity",
        shorthand: "NAE 20:4(5,8,11,14)",
        note: "N-acylethanolamine.",
        valid: true,
    },
    Example {
        group: "Headgroup diversity",
        shorthand: "AMP-FA 20:4(5,8,11,14);15OH",
        note: "AMP-linked fatty acid with a localized hydroxyl.",
        valid: true,
    },
    Example {
        group: "Headgroup diversity",
        shorthand: "Cer d18:1(4)/16:0",
        note: "Ceramide with a dihydroxy sphingoid base and N-acyl chain.",
        valid: true,
    },
    Example {
        group: "Headgroup diversity",
        shorthand: "SM d18:1(4)/16:0",
        note: "Sphingomyelin with phosphocholine headgroup.",
        valid: true,
    },
    Example {
        group: "Headgroup diversity",
        shorthand: "HexCer d18:1(4)/16:0",
        note: "Hexosylceramide with the supported hexose template.",
        valid: true,
    },
    Example {
        group: "Headgroup diversity",
        shorthand: "S1P d18:1(4)",
        note: "Sphingosine-1-phosphate.",
        valid: true,
    },
    Example {
        group: "Intentionally unsupported shorthand",
        shorthand: "PC 34:1",
        note: "Species-level sum composition does not determine two explicit chains, so no structure is invented.",
        valid: false,
    },
    Example {
        group: "Intentionally unsupported shorthand",
        shorthand: "FA 18:1(9);O2",
        note: "A generic oxygen count does not identify functional groups or candidate sites.",
        valid: false,
    },
    Example {
        group: "Intentionally unsupported shorthand",
        shorthand: "FA 18:0;1OMe",
        note: "An acyl-chain C1 substituent changes the linkage rather than modifying a chain carbon.",
        valid: false,
    },
    Example {
        group: "Intentionally unsupported shorthand",
        shorthand: "Gal-Glc-Cer d18:1(4)/16:0",
        note: "The carbohydrate sequence does not specify glycosidic linkage positions.",
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
    let canonical = generated.as_deref().and_then(canonicalize_cxsmiles);
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
        .and_then(canonicalize_cxsmiles)
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
  .notes { min-width: 260px; }
</style>
</head>
<body>
<h1>lipid_notation demo</h1>
<p>This gallery is generated by <code>tests/demo.rs</code>. Regenerate it with
<code>cargo test --test demo -- --ignored</code>.</p>
<p>Each valid shorthand is converted with <code>name2smiles</code>, canonicalized
with <code>canonicalize_cxsmiles</code>, and read back with <code>smiles2name</code>.
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
            out.push_str("<thead><tr><th>Test shorthand</th><th><code>name2smiles</code> output</th><th>Canonicalized CXSMILES</th><th><code>smiles2name</code> result</th><th>Direct CXSMILES depiction</th><th><code>smiles_for_depiction</code> depiction</th><th>Notes on correctness</th></tr></thead><tbody>\n");
        }

        let row = evaluate(example);
        validate(example, &row);
        out.push_str(&format!(
            "<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td class=\"depiction\">{}</td><td class=\"depiction\">{}</td><td class=\"notes\">{}</td></tr>\n",
            html_code(Some(example.shorthand)),
            html_code(row.generated.as_deref()),
            html_code(row.canonical.as_deref()),
            html_code(row.recovered.as_deref()),
            cdk_depiction(row.direct_depiction.as_deref()),
            cdk_depiction(row.prepared_depiction.as_deref()),
            example.note,
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
