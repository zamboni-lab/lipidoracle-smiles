use lipid_notation::{canonicalize, name2smiles, smiles2name, smiles_for_depiction};

fn main() {
    let names: &[&str] = &[
        // Fully determined (previous blog start)
        "FA 18:1(9Z)",
        // unlocalized double bond
        "FA 18:1",
        "FA 18:2",
        "FA 20:4",
        "PC 16:0/18:2",
        // modification position
        "FA 18:0;OH",
        "DG 18:0/18:2;OH(5)",
        "DG 18:0/18:2(9,12);OH",
        "DG 18:1(9);OH(5)_18:1;OH",
        // sn assignment
        "PC 16:0_18:1(9)",
        "PC 16:0_18:1",
        "DG 16:0_18:1(9)",
        "DG 16:0_18:1",
        "PC 16:0/18:1(9)",
        // consensus
        "FA 18:1(9) [DB sn1: \u{0394}9 92%]",
        "FA 18:2(9,12) [DB sn1: \u{0394}9 100%, \u{0394}12 88%]",
        "FA 20:4;OH [OH sn1: 11 50%, 13 50%]",
        "PC 16:0_18:1(9) [DB sn2: \u{0394}9 100%]",
        // demo examples
        "FA 20:4(5,8,12,14);OH",
        "PC O-16:1_18:2;OH",
        "FA 20:3(5,8,11);(OH)2",
        "FA 20:2(5,8);[11-15cy5;13OH];OH",
        "PC 16:1_18:2;9Ep;OH",
        "TG 18:1(9);5OH_18:2;9Ep_18:1",
        "Cer d18:1(4)/16:1;OH",
        // unsupported
        "PC 34:1",
        "PC 34:2",
        "TG 18:0_18:1_18:2",
    ];
    for n in names {
        let g = name2smiles(n);
        let c = g.as_deref().and_then(canonicalize);
        let r = c.as_deref().and_then(smiles2name);
        let d = g.as_deref().map(smiles_for_depiction);
        println!("=== NAME: {n}");
        println!("  generated: {}", g.as_deref().unwrap_or("<None>"));
        println!("  canonical: {}", c.as_deref().unwrap_or("<None>"));
        println!("  recovered: {}", r.as_deref().unwrap_or("<None>"));
        println!("  depiction: {}", d.as_deref().unwrap_or("<None>"));
    }
}
