use lipid_notation::{canonicalize, name2smiles, smiles2name};

/// Every (name, SMILES) pair exactly as it appears in blog/v2_followup.md.
/// 1) name2smiles must reproduce the blog SMILES byte-for-byte.
/// 2) smiles2name must read it back to some name.
/// 3) that name must regenerate a canonically equivalent structure.
static CASES: &[(&str, &str)] = &[
    (r"FA 18:1(9Z)", r"OC(=O)CCCCCCC/C=C\CCCCCCCC"),
    (r"FA 18:1", r"OC(=O)CC=CC |Sg:n:3:a:ht,Sg:n:6:b:ht| constrain(a+b=15)"),
    (r"FA 18:2", r"OC(=O)CC=CC=CC |Sg:n:3:a:ht,Sg:n:5:b:ht,Sg:n:8:c:ht| constrain(a+b+c=14)"),
    (r"FA 20:4", r"OC(=O)CC=CC=CC=CC=CC |Sg:n:3:a:ht,Sg:n:5:b:ht,Sg:n:7:c:ht,Sg:n:9:d:ht,Sg:n:12:e:ht| constrain(a+b+c+d+e=14)"),
    (r"PC 16:0/18:2", r"C(COP(=O)([O-])OCC[N+](C)(C)C)(OC(=O)CC=CC=CC)COC(=O)CCCCCCCCCCCCCCC |Sg:n:16:a:ht,Sg:n:18:b:ht,Sg:n:21:c:ht| constrain(a+b+c=14)"),
    (r"FA 18:0;OH", r"OC(=O)CCCCCCCCCCCCCCCCC.*O |m:20:3.4.5.6.7.8.9.10.11.12.13.14.15.16.17.18.19|"),
    (r"PC 16:0_18:1(9)", r"C(COP(=O)([O-])OCC[N+](C)(C)C)(OC(=O)CCCCCCCC=CCCCCCCCC)COC(=O)CCCCCCCCCCCCCCC |$;;;;;;;;;;;;;sn2;;;;;;;;;;;;;;;;;;;;;sn1$| swappable(sn1,sn2)"),
    (r"PC 16:0_18:1", r"C(COP(=O)([O-])OCC[N+](C)(C)C)(OC(=O)CC=CC)COC(=O)CCCCCCCCCCCCCCC |$;;;;;;;;;;;;;sn2;;;;;;;;sn1$,Sg:n:16:a:ht,Sg:n:19:b:ht| swappable(sn1,sn2);constrain(a+b=15)"),
    (r"DG 18:1(9);OH(5)_18:1;OH", r"C(CO)(OC(=O)CC=CC.*O)COC(=O)CCCC(O)CCCC=CCCCCCCCC |$;;;sn2;;;;;;;;;;sn1$,Sg:n:6:a:ht,Sg:n:9:b:ht,m:10:6.7.8.9| swappable(sn1,sn2);constrain(a+b=15)"),
    (r"FA 18:1(9) [DB sn1: Δ9 92%]", r"OC(=O)CCCCCCCC=CCCCCCCCC |$;sn1$| dbPos(sn1:9@92)"),
    (r"FA 20:4;OH [OH sn1: 11 50%, 13 50%]", r"OC(=O)CC=CC=CC=CC=CC.*O |$;sn1;;;;;;;;;;;;OH1$,Sg:n:3:a:ht,Sg:n:5:b:ht,Sg:n:7:c:ht,Sg:n:9:d:ht,Sg:n:12:e:ht,m:13:3.4.5.6.7.8.9.10.11.12| constrain(a+b+c+d+e=14);mPos(OH1:11OH@50,13OH@50)"),
    (r"TG 18:0_18:1_18:2", r"C(COC(=O)CC=CC=CC)(OC(=O)CC=CC)COC(=O)CCCCCCCCCCCCCCCCC |$;;sn3;;;;;;;;;sn2;;;;;;;;sn1$,Sg:n:5:a:ht,Sg:n:7:b:ht,Sg:n:10:c:ht,Sg:n:14:d:ht,Sg:n:17:e:ht| swappable(sn1,sn2,sn3);constrain(a+b+c=14);constrain(d+e=15)"),
    (r"TG 18:1(9);5OH_18:2;9Ep_18:1", r"C(COC(=O)CC=CC)(OC(=O)CCCCCCCC%10OC%10CC=CC=CC)COC(=O)CCCC(O)CCCC=CCCCCCCCC |$;;sn3;;;;;;;sn2;;;;;;;;;;;;;;;;;;;;sn1$,Sg:n:5:a:ht,Sg:n:8:b:ht,Sg:n:22:c:ht,Sg:n:24:d:ht,Sg:n:27:e:ht| swappable(sn1,sn2,sn3);constrain(a+b=15);constrain(c+d+e=5)"),
    (r"PC O-16:1_18:2;OH", r"C(COP(=O)([O-])OCC[N+](C)(C)C)(OC(=O)CC=CC=CC.*O)COCCC=CC |$;;;;;;;;;;;;;sn2;;;;;;;;;;;;sn1$,Sg:n:16:a:ht,Sg:n:18:b:ht,Sg:n:21:c:ht,Sg:n:27:d:ht,Sg:n:30:e:ht,m:22:16.17.18.19.20.21| swappable(sn1,sn2);constrain(a+b+c=14);constrain(d+e=13)"),
    (r"Cer d18:1(4)/16:1;OH", r"C(CO)(NC(=O)CC=CC.*O)C(O)C=CCCCCCCCCCCCCC |Sg:n:6:a:ht,Sg:n:9:b:ht,m:10:6.7.8.9| constrain(a+b=13)"),
];

fn main() {
    let mut pass = 0; let mut fail = 0;
    for (name, smi) in CASES {
        let gen = name2smiles(name);
        let back = smiles2name(smi);
        let exact = gen.as_deref() == Some(*smi);
        let readback = back.is_some();
        let regen_eq = match (&back, smi) {
            (Some(n), _) => canonicalize(&name2smiles(n).unwrap()) == canonicalize(smi),
            (None, _) => false,
        };
        if exact && readback && regen_eq { pass += 1; } else {
            fail += 1;
            println!("FAIL {name}: exact={exact} readback={readback} regen_eq={regen_eq} back={back:?}");
        }
    }
    println!("RESULT: {pass} passed, {fail} failed of {}", CASES.len());
    if fail > 0 { std::process::exit(1); }
}
