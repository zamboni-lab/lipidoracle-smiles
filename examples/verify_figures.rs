use lipid_notation::{name2smiles, smiles_for_depiction};
static FIGURES: &[(&str, &str, &str)] = &[
    ("fig_fa181z", r"FA 18:1(9Z)", r"C(CCCC)CCC/C=C\CCCCCCCC(=O)O"),
    ("fig_fa181", r"FA 18:1", r"CC=CCC(=O)O |Sg:n:3:a:ht,Sg:n:0:b:ht| constrain(a+b=15)"),
    ("fig_fa182", r"FA 18:2", r"C(O)(CC=CC=CC)=O |Sg:n:2:a:ht,Sg:n:4:b:ht,Sg:n:7:c:ht| constrain(a+b+c=14)"),
    ("fig_fa204", r"FA 20:4", r"C(C=CC=CC=CCC(O)=O)=CC |Sg:n:7:a:ht,Sg:n:5:b:ht,Sg:n:3:c:ht,Sg:n:1:d:ht,Sg:n:12:e:ht| constrain(a+b+c+d+e=14)"),
    ("fig_pc160_182", r"PC 16:0/18:2", r"O=P(OCC(OC(=O)CC=CC=CC)COC(CCCCCCCCCCCCCCC)=O)(OCC[N+](C)(C)C)[O-] |Sg:n:8:a:ht,Sg:n:10:b:ht,Sg:n:13:c:ht| constrain(a+b+c=14)"),
    ("fig_pc160_181_sn", r"PC 16:0_18:1(9)", r"C(CCCCCCCC)=CCCCCCCCC(OC(COC(CCCCCCCCCCCCCCC)=O)COP(=O)([O-])OCC[N+](C)(C)C)=O |$;;;;;;;;;;;;;;;;;;sn2;;;sn1$| swappable(sn1,sn2)"),
    ("fig_pc160_181", r"PC 16:0_18:1", r"O(C(=O)CCCCCCCCCCCCCCC)CC(OC(=O)CC=CC)COP([O-])(=O)OCC[N+](C)(C)C |$sn1;;;;;;;;;;;;;;;;;;;;sn2$,Sg:n:23:a:ht,Sg:n:26:b:ht| swappable(sn1,sn2);constrain(a+b=15)"),
    ("fig_fa180_oh", r"FA 18:0;OH", r"C(CCCCCCCCCCCCCCCCC)(=O)O.O[*] |m:21:1.2|"),
    ("fig_dg_combined", r"DG 18:1(9);OH(5)_18:1;OH", r"C(OC(COC(CCCC(CCCC=CCCCCCCCC)O)=O)CO)(CC=CC)=O.O[*] |$;sn2;;;sn1$,Sg:n:27:a:ht,Sg:n:30:b:ht,m:33:27.28| swappable(sn1,sn2);constrain(a+b=15)"),
    ("fig_fa181_conf", r"FA 18:1(9) [DB sn1: Δ9 92%]", r"C(CCCC)CCCC=CCCCCCCCC(=O)O |$;;;;;;;;;;;;;;;;;sn1$| dbPos(sn1:9@92)"),
    ("fig_fa204_oh_conf", r"FA 20:4;OH [OH sn1: 11 50%, 13 50%]", r"C(C=CC=CC=CCC(O)=O)=CC.O[*] |$;;;;;;;;sn1;;;;;;OH1$,Sg:n:7:a:ht,Sg:n:5:b:ht,Sg:n:3:c:ht,Sg:n:1:d:ht,Sg:n:12:e:ht,m:14:7.6| constrain(a+b+c+d+e=14);mPos(OH1:11OH@50,13OH@50)"),
    ("fig_tg_3chain", r"TG 18:0_18:1_18:2", r"CC=CCC(=O)OC(COC(=O)CC=CC=CC)COC(CCCCCCCCCCCCCCCCC)=O |$;;;;;;sn2;;;sn3;;;;;;;;;;sn1$,Sg:n:12:a:ht,Sg:n:14:b:ht,Sg:n:17:c:ht,Sg:n:3:d:ht,Sg:n:0:e:ht| swappable(sn1,sn2,sn3);constrain(a+b+c=14);constrain(d+e=15)"),
    ("fig_demo_tg_epoxide", r"TG 18:1(9);5OH_18:2;9Ep_18:1", r"C(=O)(CCCCCCCC1C(CC=CC=CC)O1)OC(COC(CC=CC)=O)COC(CCCC(CCCC=CCCCCCCCC)O)=O |$;;;;;;;;;;;;;;;;;;sn2;;;sn3;;;;;;;;sn1$,Sg:n:23:a:ht,Sg:n:26:b:ht,Sg:n:11:c:ht,Sg:n:13:d:ht,Sg:n:16:e:ht| swappable(sn1,sn2,sn3);constrain(a+b=15);constrain(c+d+e=5)"),
    ("fig_demo_pc_ether", r"PC O-16:1_18:2;OH", r"O[*].C(C)=CCCOCC(OC(CC=CC=CC)=O)COP(=O)([O-])OCC[N+](C)(C)C |$;;;;;;;sn1;;;sn2$,Sg:n:12:a:ht,Sg:n:14:b:ht,Sg:n:17:c:ht,Sg:n:5:d:ht,Sg:n:3:e:ht,m:1:12.13| swappable(sn1,sn2);constrain(a+b+c=14);constrain(d+e=13)"),
    ("fig_demo_cer", r"Cer d18:1(4)/16:1;OH", r"C(NC(C(O)C=CCCCCCCCCCCCCC)CO)(CC=CC)=O.O[*] |Sg:n:22:a:ht,Sg:n:25:b:ht,m:28:22.23| constrain(a+b=13)"),
];
fn main(){ let mut pass=0; let mut fail=0;
  for (f, name, expected) in FIGURES {
    let got = name2smiles(name).map(|g| smiles_for_depiction(&g));
    let ok = got.as_deref() == Some(*expected);
    if ok { pass+=1; } else { fail+=1; println!("FAIL {f}: expected {:?} got {:?}", expected, got); }
  }
  println!("RESULT: {pass} passed, {fail} failed of {}", FIGURES.len());
  if fail>0 { std::process::exit(1); }
}
