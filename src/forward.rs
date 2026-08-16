//! Forward conversion from lipid notation to SMILES/CXSMILES.
//!
//! `generate_smiles` takes a canonical Shorthand2020 name —
//! e.g. `"TG 18:0/18:0/18:1(9);5oxo"`,
//! `"PC 16:0/20:4(5,8,11,14);15OH"`, or the documented project extension
//! `"AMP-FA 20:4(5,8,11,14);15OH"` — and builds a flat (non-stereo)
//! SMILES for that isomer. Double-bond/hydroxyl/ketone/extra-carboxyl
//! positions are carbon indices counted from C1 (the carboxyl carbon for
//! acyl chains, C1 of the sphingoid base otherwise) — the same Δ-numbering
//! LipidOracle's EAD engines use in their fragment labels.
//!
//! ## Shorthand vs explicit chains
//!
//! Multi-chain lipid **shorthand** (e.g., `"PC 19:2"` or `"TG 54:3"`) — total
//! composition without explicit chain breakdown — is rejected and returns `None`.
//! Such shorthand is ambiguous: `"PC 19:2"` could represent many chain combinations
//! (16:0/3:2, 15:1/4:1, etc.), and there's no canonical single-chain representation.
//! Supply **explicit chains** instead: e.g., `"PC 16:0/18:2"` or `"PC 16:0_18:2"`.
//! Single-chain lipids (FA, CE, LPC, LPE, etc.) can use shorthand—e.g., `"FA 18:1"`.
//!
//! ## Separator semantics
//!
//! `/` between explicit chains means sn-position is known: chains are read
//! positionally (sn1, sn2, sn3, ...) and rendered as one connected
//! molecule. `_` means regiochemistry is *not* known (the classic
//! bulk/species-level lipid shorthand, and also what LipidOracle's EAD
//! engines always use for DG/TG/CL regardless of how well the individual
//! chains are localized) — in that case the chains are still written into
//! the backbone in name order, but the atom each one hangs from is labelled `snN`
//! and a `swappable(...)` token after
//! the closing pipe says their assignment is one arbitrary choice among
//! several, e.g. `PC 16:0_18:1(9)` becomes
//! `C(COP(=O)([O-])OCC[N+](C)(C)C)(OC(=O)CCCCCCCC=CCCCCCCCC)COC(=O)CCCCCCCCCCCCCCC |$;;;;;;;;;;;;;sn2;;;;;;;;;;;;;;;;;;;;;sn1$| swappable(sn1,sn2)`.
//!
//! ## Unlocalized double bonds
//!
//! Every chain builder is Sg:-aware: known double bonds/modifications are
//! rendered literally (with geometry where declared), and any remaining
//! *unlocalized* double bonds are represented with a CDK-style `Sg:n:`
//! flexible-run marker plus an `x+y=N` size constraint, never with a
//! guessed literal position. See `build_chain_segment` and
//! `CxBuilder` for how a chain's local Sg:/m: blocks get offset and
//! merged into the surrounding headgroup template's global CXSMILES
//! suffix. A double bond whose geometry is not explicit is left as a plain
//! `C=C`, which is exactly what unspecified geometry means in SMILES;
//! explicit `Z`/`E` geometry remains encoded with
//! `/` and `\\`.
//!
//! ## Supported modifications
//!
//! Every functional group in Table 1A of Liebisch et al. 2020 (PMC7707175)
//! that substitutes a chain carbon, and every `cyX` ring in Table 1B. See
//! [`SUBSTITUENTS`] for the branch each group renders as, and [`Ring`] for
//! the ring spellings.
//!
//! Positions go in front of the abbreviation and are comma-separated
//! (`;11OH,15OH;9oxo`); a group named without one (`;OH`) is declared
//! present but unlocalized and becomes an `m:` block; the parenthesized
//! form carries a count or protects an abbreviation containing digits
//! (`;(OH)2`, `;(NO2)`). Compatibility aliases such as `;OH(3,5)` and
//! `;ep(5)` are also accepted.
//!
//! `None` is returned when:
//!
//! * the headgroup has no template here;
//! * a chain's oxygen count is declared generically (`;O`/`;O2`) with no
//!   group breakdown at all — unlike double bonds, no placeholder
//!   convention was asked for there;
//! * the group is one of [`UNRENDERABLE`], or sits on C1 of an acyl chain,
//!   where the shorthand is naming a linkage rather than a substituent;
//! * a chain's unlocalized double bonds have less room left than the `Sg:`
//!   scaffold needs (see `build_chain_segment`).
//!
//! Table 1C's carbohydrates are headgroups rather than chain
//! modifications — `Gal-Glc-Cer` names the sugar *sequence* but never the
//! glycosidic linkage positions, so there is no honest single structure to
//! emit and they are not handled here.

use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Geometry {
    Cis,
    Trans,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ChainPrefix {
    Acyl,
    EtherAlkyl,
    EtherAlkenyl,
    SphingoidD,
    SphingoidT,
}

/// One double-bond position: the Δ-carbon it starts at, its geometry (if
/// given), and whether the position itself is a guessed placeholder
/// (declared count exceeded what was actually localized).
#[derive(Debug, Clone, Copy)]
struct DbPos {
    pos: u32,
    geom: Option<Geometry>,
    placeholder: bool,
}

/// Table 1A of Liebisch et al. 2020 (PMC7707175), in the paper's own
/// IUPAC-hierarchy order: every functional group that substitutes a single
/// chain carbon, paired with the SMILES branch that carbon then carries.
/// Carbon `k` is emitted as `C(<branch>)`, so `oxo`'s branch is `=O` and
/// its carbon becomes the `C(=O)` a ketone needs.
///
/// Matching is on the whole token, never a prefix, so `OOH` and `OH`, or
/// `S` and `SH`, cannot be confused for one another.
pub(crate) const SUBSTITUENTS: &[(&str, &str)] = &[
    ("Et", "CC"),
    ("Me", "C"),
    ("Br", "Br"),
    ("Cl", "Cl"),
    ("F", "F"),
    ("I", "I"),
    ("NO2", "[N+](=O)[O-]"),
    ("OMe", "OC"),
    ("NH2", "N"),
    ("OOH", "OO"),
    ("SH", "S"),
    ("OH", "O"),
    ("oxo", "=O"),
    ("CN", "C#N"),
    ("P", "OP(=O)(O)O"),
    ("S", "OS(=O)(=O)O"),
    ("COOH", "C(=O)O"),
];

/// Table 1A entries this generator refuses rather than guesses.
///
/// `oxy` (alkoxy) and `OO` (peroxy) name a linkage whose other end the
/// shorthand never states — an ether to *what*, a peroxide bridging *which*
/// second carbon — and `G`/`T` conjugate glycine or taurine to the carboxyl
/// terminus, which is a change of headgroup rather than a substituent on a
/// chain carbon. Emitting any of them would mean inventing the missing half,
/// which is the one thing this crate exists not to do.
const UNRENDERABLE: &[&str] = &["oxy", "OO", "G", "T"];

/// One functional group on a chain: its Δ-position (`0` meaning "declared
/// present, position not determined") and its `SUBSTITUENTS` branch.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Mod {
    pos: u32,
    branch: &'static str,
}

/// A ring closed across a span of chain carbons: the `cyX` carbocycles of
/// Table 1B (`FA 19:1;[11-13cy3:0]`, lactobacillic acid), or an `Ep` epoxide,
/// which is the same thing with an oxygen sitting in the ring.
///
/// `start`/`end` are the chain carbons the ring closure bonds together; for
/// an epoxide they are adjacent and `bridge` supplies the ring's third atom.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Ring {
    start: u32,
    end: u32,
    bridge: Option<&'static str>,
}

#[derive(Debug, Clone)]
struct ParsedChain {
    prefix: ChainPrefix,
    carbon: u32,
    db_pos: Vec<DbPos>,
    mods: Vec<Mod>,
    rings: Vec<Ring>,
}

/// Entry point: build semantically correct SMILES/CXSMILES for one lipid species.
///
/// For structures with variable content (unknown DBs, unknown mods, or
/// unresolved regiochemistry), returns CXSMILES with `Sg:`, `m:`, and atom
/// labels. Fully determined structures return plain SMILES. Unsupported or
/// structurally under-specified names return `None`.
///
/// Reference: https://egonw.github.io/cdk-cxsmiles/templates.html#lipids-with-two-double-bonds-somewhere-in-the-tail
pub(crate) fn generate_smiles(name: &str) -> Option<String> {
    let (canonical, tail) = crate::nomenclature::split_display_name(name);
    let built = build_lipid(&canonical, true)?;
    Some(
        match tail.as_deref().and_then(crate::consensus::parse_tail) {
            Some(entries) => attach_consensus(built, &entries),
            None => built.smiles,
        },
    )
}

/// Adds the `dbPos`/`mPos` tokens for a name's consensus tail, plus whatever
/// anchors they need.
///
/// The tokens name atom labels rather than atom indices, so this makes sure
/// the labels exist: the `$snN$` block is emitted even when the structure had
/// no other reason for one, and each `m:` stub is labelled so a token about a
/// floating group can point at that group rather than merely at its chain.
///
/// Nothing here touches the CXSMILES fields themselves — the consensus is
/// metadata and stays in the trailer, which is what it is for.
fn attach_consensus(built: Built, entries: &[crate::consensus::Consensus]) -> String {
    let (base, fields, trailer) = match built.smiles.split_once(" |") {
        Some((base, rest)) => match rest.split_once('|') {
            Some((fields, trailer)) => (base, fields.to_string(), trailer.trim().to_string()),
            None => (built.smiles.as_str(), String::new(), String::new()),
        },
        None => (built.smiles.as_str(), String::new(), String::new()),
    };

    // Label each `m:` stub `<abbr><n>`, so an unlocalized group has a name.
    let stubs = stub_labels(base, &fields, &built);
    let mut sites: Vec<(String, usize)> = built
        .sn_sites
        .iter()
        .map(|(sn, atom)| (format!("sn{sn}"), *atom))
        .collect();
    sites.extend(
        stubs
            .iter()
            .map(|(label, atom, _, _)| (label.clone(), *atom)),
    );
    sites.sort_by_key(|(_, atom)| *atom);

    let tokens = crate::consensus::tokens(entries, |entry| {
        // A token about a group the structure left floating anchors to that
        // stub; one about a position the structure already states anchors to
        // the chain.
        (!entry.is_localized())
            .then(|| {
                stubs
                    .iter()
                    .find(|(_, _, kind, sn)| *kind == entry.kind && *sn == entry.sn)
                    .map(|(label, _, _, _)| label.clone())
            })
            .flatten()
    });

    let fields = match fields.contains('$') {
        true => fields,
        false => {
            let block = label_block(&sites);
            if fields.is_empty() {
                block
            } else {
                format!("{block},{fields}")
            }
        }
    };

    let mut trailer_tokens: Vec<String> = trailer
        .split(';')
        .filter(|t| !t.trim().is_empty())
        .map(str::to_string)
        .collect();
    trailer_tokens.extend(tokens);
    format!("{base} |{fields}| {}", trailer_tokens.join(";"))
        .trim_end()
        .to_string()
}

/// `(label, atom, group kind, sn)` for every `m:` stub in a built structure.
///
/// The stub's own component says which group it is (`*O` is a hydroxyl), and
/// its candidate list says which chain, so a consensus entry can be matched to
/// the stub it is about.
fn stub_labels(base: &str, fields: &str, built: &Built) -> Vec<(String, usize, String, usize)> {
    let components: Vec<&str> = base.split('.').collect();
    let mut out = Vec::new();
    let mut counts: HashMap<String, usize> = HashMap::new();

    for field in fields.split(',') {
        let Some(rest) = field.strip_prefix("m:") else {
            continue;
        };
        let Some((atom, candidates)) = rest.split_once(':') else {
            continue;
        };
        let Ok(atom) = atom.parse::<usize>() else {
            continue;
        };
        let sites: Vec<usize> = candidates
            .split('.')
            .filter_map(|s| s.parse().ok())
            .collect();
        let Some((sn, _)) = built
            .chains
            .iter()
            .find(|(_, carbons)| sites.iter().all(|s| carbons.contains(s)))
        else {
            continue;
        };
        // The stub's branch, e.g. `*O`, names the group via SUBSTITUENTS.
        let branch = components
            .iter()
            .skip(1)
            .find_map(|c| c.strip_prefix('*'))
            .unwrap_or("");
        let Some((kind, _)) = SUBSTITUENTS.iter().find(|(_, b)| *b == branch) else {
            continue;
        };
        let n = counts.entry((*kind).to_string()).or_insert(0);
        *n += 1;
        out.push((format!("{kind}{n}"), atom, (*kind).to_string(), *sn));
    }
    out
}

/// A `$...$` block from `(label, atom)` pairs.
fn label_block(sites: &[(String, usize)]) -> String {
    let last = sites.iter().map(|(_, atom)| *atom).max().unwrap_or(0);
    let labels: Vec<&str> = (0..=last)
        .map(|i| {
            sites
                .iter()
                .find(|(_, atom)| *atom == i)
                .map(|(label, _)| label.as_str())
                .unwrap_or("")
        })
        .collect();
    format!("${}$", labels.join(";"))
}

/// Builds one lipid's SMILES plus, for every chain, the atom index of each
/// of its carbons — the mapping a UI needs to highlight the part of a
/// structure that a given MS2 fragment comes from.
///
/// Differs from [`generate_smiles`] in two depiction-driven ways:
///
/// * `_`-joined names use one representative chain assignment and set
///   [`LipidStructure::regio_resolved`] to `false`;
/// * The `Sg:` unlocalized-double-bond markers are already expanded (see
///   [`smiles_expand`]), and the atom indices account for
///   the padding atoms that expansion inserts.
///
/// Atom indices are 0-based in SMILES emission order, which is the order
/// RDKit (and every other parser) assigns when reading the string.
pub(crate) fn generate_structure(name: &str) -> Option<LipidStructure> {
    let (canonical, _) = crate::nomenclature::split_display_name(name);
    let labels = chain_label_tokens(&canonical);
    let built = build_lipid(&canonical, false)?;

    let (expanded, inserts) = expand_with_padding_inserts(&built.smiles);
    let atom_count = count_atoms(&expanded);
    let mut shift = inserts.clone();
    shift.sort_by_key(|(pos, _)| *pos);
    let insert_counts: HashMap<usize, usize> = inserts.iter().copied().collect();

    let chains = built
        .chains
        .into_iter()
        .filter(|(_, carbons)| !carbons.is_empty())
        .map(|(sn, carbons)| {
            let mut mapped = Vec::with_capacity(carbons.len());
            for atom in carbons {
                let base = shift_atom(atom, &shift);
                mapped.push(base);
                // Expansion writes each marker's padding atoms directly
                // after it, so they continue this chain's carbon run.
                for k in 1..=insert_counts.get(&atom).copied().unwrap_or(0) {
                    mapped.push(base + k);
                }
            }
            ChainAtoms {
                sn,
                label: labels.get(sn - 1).cloned().unwrap_or_default(),
                carbons: mapped,
            }
        })
        .collect();

    Some(LipidStructure {
        smiles: expanded,
        chains,
        atom_count,
        regio_resolved: !is_sn_unresolved(&canonical),
    })
}

/// Per-chain atom indices for one chain of a [`LipidStructure`].
#[derive(Debug, Clone)]
pub struct ChainAtoms {
    /// 1-based sn position, counted over the name's chain tokens (including
    /// any `0:0` placeholder), so `sn` stays aligned with the `snN` prefix
    /// the EAD engines put on their fragment labels.
    pub sn: usize,
    /// The chain's own token from the lipid name, e.g. `"18:1(9)"`.
    pub label: String,
    /// Atom index of C1, C2, ... Cn, in chain numbering order (C1 is the
    /// carboxyl carbon for acyl chains, C1 of the base for sphingoids —
    /// the same numbering the EAD `C{k}-C{k+1}` fragment labels use).
    pub carbons: Vec<usize>,
}

/// A depiction-ready structure plus the atom indices of every chain carbon.
#[derive(Debug, Clone)]
pub struct LipidStructure {
    /// Plain SMILES, already expanded — no CXSMILES suffix.
    pub smiles: String,
    pub chains: Vec<ChainAtoms>,
    pub atom_count: usize,
    /// False when the source name used `_` and the sn assignment shown is
    /// therefore one arbitrary choice among several.
    pub regio_resolved: bool,
}

/// The chain tokens of a lipid name, in sn order, `/` and `_` alike.
fn chain_label_tokens(name: &str) -> Vec<String> {
    let Some((_, rest)) = name.trim().split_once(' ') else {
        return Vec::new();
    };
    rest.trim()
        .replace('/', "_")
        .split('_')
        .map(|s| s.trim().to_string())
        .collect()
}

fn is_sn_unresolved(name: &str) -> bool {
    match name.trim().split_once(' ') {
        Some((_, rest)) => !rest.contains('/') && rest.contains('_'),
        None => false,
    }
}

/// Shared builder behind [`generate_smiles`] and
/// [`generate_structure`]. `force` is the regiochemistry mode to use
/// for multi-chain `_`-joined names; `/`-joined and single-chain names are
/// always resolved regardless.
fn build_lipid(name: &str, mark_swappable: bool) -> Option<Built> {
    let (class, rest) = match name.trim().split_once(' ') {
        Some((h, r)) => (h, r.trim()),
        None => (name.trim(), ""),
    };

    if rest.is_empty() {
        // Single lipid class with no chains - return simple SMILES
        return (class == "ST").then(|| Built::plain(st_smiles()));
    }

    let is_slash = rest.contains('/');
    let has_separator = is_slash || rest.contains('_');
    let normalized = rest.replace('/', "_");
    let chains: Vec<ParsedChain> = normalized
        .split('_')
        .map(parse_chain_token)
        .collect::<Option<Vec<_>>>()?;

    // Reject shorthand notation (no separators) when headgroup expects multiple chains
    if !has_separator && chains.len() == 1 && class_needs_multi_chain(class) {
        return None;
    }

    // Only a `_`-joined multi-chain name has an sn ambiguity to state.
    let swap = mark_swappable && !is_slash && chains.len() > 1;

    match class {
        "FA" => fa_smiles(&chains),
        "AMP-FA" | "FA-AMP" => amp_fa_smiles(&chains),
        "NAE" => nae_smiles(&chains),
        "CAR" => car_smiles(&chains),
        "CE" => ce_smiles(&chains),
        "ST" => Some(Built::plain(st_smiles())),
        "MG" => glycerolipid_smiles(&chains, 1, false),
        "DG" => glycerolipid_smiles(&chains, 2, swap),
        "TG" => glycerolipid_smiles(&chains, 3, swap),
        "CL" => cl_smiles(&chains, swap),
        "PC" | "LPC" => gpl_smiles("OP(=O)([O-])OCC[N+](C)(C)C", &chains, swap),
        "PE" | "LPE" => gpl_smiles("OP(=O)(O)OCCN", &chains, swap),
        "PS" | "LPS" => gpl_smiles("OP(=O)(O)OCC(N)C(=O)O", &chains, swap),
        "PG" | "LPG" => gpl_smiles("OP(=O)(O)OCC(O)CO", &chains, swap),
        "PI" | "LPI" => gpl_smiles("OP(=O)(O)OC1C(O)C(O)C(O)C(O)C1O", &chains, swap),
        "PA" | "LPA" => gpl_smiles("OP(=O)(O)O", &chains, swap),
        "Cer" => sphingo_dispatch(&chains, "O", true),
        "CerP" => sphingo_dispatch(&chains, "OP(=O)(O)O", true),
        "SM" => sphingo_dispatch(&chains, "OP(=O)([O-])OCC[N+](C)(C)C", true),
        "HexCer" => sphingo_dispatch(&chains, "OC1OC(CO)C(O)C(O)C1O", true),
        "IPC" => sphingo_dispatch(&chains, "OP(=O)(O)OC1C(O)C(O)C(O)C(O)C1O", true),
        "S1P" => sphingo_dispatch(&chains, "OP(=O)(O)O", false),
        "Sph" | "SB" => sphingo_dispatch(&chains, "O", false),
        _ => None,
    }
}

/// Classes whose structure inherently spans multiple acyl/alkyl chains
/// (sn1/sn2[/sn3/sn4]), so a bare sum-composition token (e.g. `"32:1"`,
/// no `"_"`/`"/"` chain split) can't be resolved to one real structure —
/// the total carbons/double bonds could come from many different chain
/// combinations. Used both to reject shorthand species names here and,
/// by the EAD engines, to reject the same shorthand as a double-bond
/// localization target in the first place (localizing a position within
/// the *sum* would treat two real chains as if they were one).
pub(crate) fn class_needs_multi_chain(class: &str) -> bool {
    matches!(
        class,
        "DG" | "TG"
            | "CL"
            | "PC"
            | "PE"
            | "PS"
            | "PG"
            | "PI"
            | "PA"
            | "Cer"
            | "CerP"
            | "SM"
            | "HexCer"
            | "IPC"
    )
}

/// An assembled structure: the SMILES text plus, per chain, the sn position
/// and the local atom indices of its carbons in C1..Cn order.
#[derive(Debug, Default, Clone)]
struct Built {
    smiles: String,
    chains: Vec<(usize, Vec<usize>)>,
    /// `(sn, atom)` of each chain's linking atom, for anchoring trailer
    /// tokens that name a chain.
    sn_sites: Vec<(usize, usize)>,
}

impl Built {
    /// A structure with no chain mapping (fixed templates like sterol).
    fn plain(smiles: String) -> Self {
        Built {
            smiles,
            chains: Vec::new(),
            sn_sites: Vec::new(),
        }
    }
}

/// A chain fragment with CXSMILES ambiguity annotations, using local atom
/// indices (0-based, counting every
/// element/bracket-atom/`*` token from this fragment's own first atom)
/// and LOCAL variable letters (always starting fresh at `a`). Embed into
/// a larger assembled SMILES via `CxBuilder::push_fragment`, which
/// offsets positions and renames variables so they don't collide with
/// whatever came before in the assembly.
#[derive(Debug, Default, Clone)]
struct CxFragment {
    smiles: String,
    sg_blocks: Vec<String>,
    constraint: Option<String>,
    m_blocks: Vec<String>,
    /// Local atom index of each chain carbon actually emitted, in C1..Cn
    /// order. Shorter than the chain's carbon count when `Sg:` markers
    /// stand in for a variable-length run — the atoms that
    /// `smiles_expand` inserts fill the rest of the run in
    /// place, so `generate_structure` splices them back in after
    /// expansion. Empty for non-chain fragments (bare `O` slots).
    carbon_atoms: Vec<usize>,
}

/// Builds one chain's fragment (C1..Cn), dispatching on its linkage type.
/// `EtherAlkenyl` gets its mandatory (never placeholder) vinyl-ether
/// C1=C2 folded in before any of the chain's own declared double bonds.
/// Sphingoid bases aren't chain fragments in this sense — see
/// `sphingoid_smiles`, which uses `build_chain_segment` directly
/// with a `start` of 3 (after the fixed C1/C2 positions).
fn build_chain_fragment(chain: &ParsedChain) -> Option<CxFragment> {
    match chain.prefix {
        ChainPrefix::Acyl => build_chain_segment(
            1,
            chain.carbon,
            &chain.db_pos,
            &chain.mods,
            &chain.rings,
            true,
        ),
        ChainPrefix::EtherAlkyl => build_chain_segment(
            1,
            chain.carbon,
            &chain.db_pos,
            &chain.mods,
            &chain.rings,
            false,
        ),
        ChainPrefix::EtherAlkenyl => {
            let mut db = vec![DbPos {
                pos: 1,
                geom: None,
                placeholder: false,
            }];
            db.extend(chain.db_pos.iter().copied());
            build_chain_segment(1, chain.carbon, &db, &chain.mods, &chain.rings, false)
        }
        ChainPrefix::SphingoidD | ChainPrefix::SphingoidT => None,
    }
}

/// Builds a chain fragment spanning positions `start..=carbon`: known
/// double bonds/modifications (`pos != 0`, `placeholder == false`) are
/// rendered literally, with geometry where declared. Any remaining
/// *unlocalized* double bonds (however many `db_pos` entries have
/// `placeholder == true` — their own guessed `.pos` is discarded, only
/// the count matters) are represented as an `Sg:n:` flexible run
/// spanning everything after the last known feature, plus a size
/// constraint, rather than a guessed literal position: each Sg-marked
/// atom is itself the first unit of its own variable's count, so
/// `smiles_expand` only needs to insert `value - 1` more
/// to reach any one valid total length. Modifications with no position
/// at all (`pos == 0`) become extra dot-joined fragments plus `m:`
/// blocks, appended after the main chain.
fn build_chain_segment(
    start: u32,
    carbon: u32,
    db_pos: &[DbPos],
    mods: &[Mod],
    rings: &[Ring],
    carbonyl_c1: bool,
) -> Option<CxFragment> {
    if carbon < start {
        return Some(CxFragment::default());
    }

    let localized: Vec<DbPos> = db_pos.iter().copied().filter(|d| !d.placeholder).collect();
    let unlocalized_count = db_pos.len() - localized.len();
    let known_mods: Vec<Mod> = mods.iter().copied().filter(|m| m.pos != 0).collect();

    // On an acyl chain C1 is the ester/amide carbonyl. A group written
    // there (`FA 18:0;1OMe`, the methyl ester) is describing the linkage,
    // not a substituent on a chain carbon, and this builder has no way to
    // change the headgroup it is being attached to. Refusing beats emitting
    // the unmodified chain and losing the group without a word.
    if carbonyl_c1 && start == 1 && known_mods.iter().any(|m| m.pos == 1) {
        return None;
    }

    // A localized modification or ring pins a chain carbon and splits the
    // chain into a proximal and a distal segment. The flexible run below only
    // ever covers what follows the last pinned feature, so every unlocalized
    // double bond would be committed to the distal segment — a claim the name
    // never made. Shorthand has no syntax for how the bonds distribute across
    // the split (`18:2;9OH` states a total and nothing else), so there is no
    // name to fall back to and no honest string to emit. Refusing beats
    // silently narrowing the name.
    //
    // A pin close enough to the start leaves no room for a double bond before
    // it, so only one distribution is consistent and the string is exact.
    if unlocalized_count > 0 {
        let first_pin = known_mods
            .iter()
            .map(|m| m.pos)
            .chain(rings.iter().map(|r| r.start))
            .min();
        // C1 of an acyl chain is the carbonyl and can never carry a double
        // bond, so the first bond-capable carbon sits one further along.
        let first_free = if carbonyl_c1 { start + 1 } else { start };
        if first_pin.is_some_and(|pin| pin > first_free + 1) {
            return None;
        }
    }

    let mut smiles;
    let mut sg_blocks = Vec::new();
    let mut constraint = None;
    let mut carbon_atoms;

    if unlocalized_count == 0 {
        (smiles, carbon_atoms) = build_chain_range(
            carbon,
            &localized,
            &known_mods,
            rings,
            carbonyl_c1,
            start,
            carbon,
        );
    } else {
        // The literal prefix has to run far enough to carry every feature
        // whose position *is* known; the flexible run covers only what is
        // left after the last of them.
        let mut prefix_len = start;
        for d in &localized {
            prefix_len = prefix_len.max(d.pos + 1);
        }
        for m in &known_mods {
            prefix_len = prefix_len.max(m.pos);
        }
        for r in rings {
            prefix_len = prefix_len.max(r.end);
        }
        prefix_len = prefix_len.min(carbon);

        (smiles, carbon_atoms) = build_chain_range(
            carbon,
            &localized,
            &known_mods,
            rings,
            carbonyl_c1,
            start,
            prefix_len,
        );

        let remaining_length = carbon - prefix_len;
        // The fixed scaffold contains `2 + 2n` carbons. Expansion only adds
        // atoms, so shorter tails cannot represent the declared chain length.
        if remaining_length < 2 * unlocalized_count as u32 + 2 {
            return None;
        }
        let var_names = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j'];
        let num_markers = (unlocalized_count + 1).min(var_names.len());
        let prefix_atom_count = count_atoms(&smiles);

        smiles.push('C');
        for _ in 0..unlocalized_count.saturating_sub(1) {
            smiles.push_str("C=C");
        }
        // Keep a terminal saturated carbon after the last uncertain double
        // bond.  Besides avoiding a chain that ends in C=C, this gives the
        // final Sg variable a real terminal atom to sit on.
        smiles.push_str("C=CC");
        // Every trailing atom is itself a chain carbon, continuing the
        // C1..Cn run straight on from the prefix.
        for i in 0..(2 + 2 * unlocalized_count) {
            carbon_atoms.push(prefix_atom_count + i);
        }

        for (i, &var) in var_names.iter().enumerate().take(num_markers) {
            let atom = if i + 1 == num_markers {
                // The last variable is anchored on the terminal carbon of
                // the C=CC suffix; the other variables begin at every other
                // carbon in the preceding flexible run.
                prefix_atom_count + 1 + 2 * unlocalized_count
            } else {
                prefix_atom_count + i * 2
            };
            sg_blocks.push(format!("Sg:n:{}:{}:ht", atom, var));
        }
        // Fixed scaffold atoms are excluded from the variable-length total.
        let constraint_sum = remaining_length.saturating_sub(unlocalized_count as u32 + 1);
        constraint = Some(format!(
            "{}={}",
            var_names[..num_markers]
                .iter()
                .map(|v| v.to_string())
                .collect::<Vec<_>>()
                .join("+"),
            constraint_sum
        ));
    }

    let mut m_blocks = Vec::new();
    let unlocalized_mods: Vec<Mod> = mods.iter().copied().filter(|m| m.pos == 0).collect();
    if !unlocalized_mods.is_empty() {
        // A CXSMILES position-variation block is `m:<floating atom>:<the
        // atoms it may be attached to>`, so the candidate list has to be
        // this chain's own carbons — the modification is declared present
        // but unlocalized, meaning *any* of them. C1 is excluded on acyl
        // chains: it's the ester/amide carbonyl, not a substitutable
        // methylene. (On ether/sphingoid ranges the first emitted carbon
        // is a legitimate site, so nothing is skipped there.)
        let skip = usize::from(carbonyl_c1 && start == 1);
        let sites = carbon_atoms
            .iter()
            .skip(skip)
            .map(|a| a.to_string())
            .collect::<Vec<_>>()
            .join(".");
        if sites.is_empty() {
            return None; // nowhere on this chain to put the modification
        }
        let mut atom_count = count_atoms(&smiles);
        // Position-variation bonds target a singly bonded wildcard. The
        // component is the wildcard plus the substituent branch (`*O`,
        // `*=O`, `*C(=O)O`, and so on).
        for m in &unlocalized_mods {
            let component = format!("*{}", m.branch);
            smiles.push('.');
            smiles.push_str(&component);
            m_blocks.push(format!("m:{atom_count}:{sites}"));
            atom_count += count_atoms(&component);
        }
    }

    Some(CxFragment {
        smiles,
        sg_blocks,
        constraint,
        m_blocks,
        carbon_atoms,
    })
}

/// A glycerol-ester/ether slot: `"O" + fragment` for a real chain, or a
/// bare free hydroxyl `"O"` when the position is empty/absent.
fn build_chain_slot(chain: Option<&ParsedChain>) -> Option<CxFragment> {
    match chain {
        None => Some(CxFragment {
            smiles: "O".to_string(),
            ..Default::default()
        }),
        Some(c) if c.carbon == 0 => Some(CxFragment {
            smiles: "O".to_string(),
            ..Default::default()
        }),
        Some(c) => {
            let frag = build_chain_fragment(c)?;
            // Prepending "O" adds one atom ahead of everything in `frag`,
            // so its own Sg:/m: positions (local to frag.smiles alone)
            // must shift by 1 to stay correct in the combined smiles —
            // otherwise they land one atom too early (e.g. inside the
            // C1 carbonyl's own `(=O)` branch instead of after it).
            let no_renaming = HashMap::new();
            let sg_blocks = frag
                .sg_blocks
                .iter()
                .map(|b| offset_sg_block(b, 1, &no_renaming))
                .collect();
            let m_blocks = frag.m_blocks.iter().map(|b| offset_m_block(b, 1)).collect();
            Some(CxFragment {
                smiles: format!("O{}", frag.smiles),
                sg_blocks,
                constraint: frag.constraint,
                m_blocks,
                carbon_atoms: frag.carbon_atoms.iter().map(|a| a + 1).collect(),
            })
        }
    }
}

/// Incrementally assembles a multi-part CXSMILES body (fixed literal
/// text interleaved with chain fragments), offsetting each fragment's
/// local `Sg:`/`m:` atom indices by the atom count of everything emitted
/// before it, and renaming each fragment's own local variable letters
/// (which always start fresh at `a`) so they don't collide with an
/// earlier fragment's.
#[derive(Default)]
struct CxBuilder {
    smiles: String,
    sg_blocks: Vec<String>,
    m_blocks: Vec<String>,
    constraints: Vec<String>,
    /// `(sn, atom index)` of the linking atom of each chain pushed via
    /// `push_chain` — the ester/ether oxygen a chain hangs from. Rendered
    /// as `$...;sn1;...$` atom labels when `swappable` is set, which is what
    /// gives the `swappable(...)` token something in the string to name.
    sn_sites: Vec<(usize, usize)>,
    /// Whether the sn assignment shown is one arbitrary choice among
    /// several, i.e. the name joined its chains with `_`.
    swappable: bool,
    atom_offset: usize,
    var_offset: usize,
    /// How many ring-closure labels earlier fragments have already claimed.
    /// Each chain numbers its own rings from `%10`; this shifts them so two
    /// chains in the same molecule cannot close each other's rings.
    ring_offset: usize,
    /// `(sn, global atom indices of C1..Cn)` for each chain pushed via
    /// `push_chain`. Recorded explicitly rather than by push order, since
    /// several builders emit their chains out of sn order (`gpl_smiles`
    /// writes sn2 before sn1).
    chains: Vec<(usize, Vec<usize>)>,
}

impl CxBuilder {
    /// Appends fixed, unambiguous literal SMILES text (headgroup
    /// backbone pieces, linking atoms, ring closures, ...).
    fn push_fixed(&mut self, text: &str) {
        self.atom_offset += count_atoms(text);
        self.smiles.push_str(text);
    }

    /// Appends a chain fragment, offsetting/renaming its local `Sg:`/`m:`
    /// blocks and constraint into this assembly's global numbering.
    fn push_fragment(&mut self, frag: &CxFragment) {
        let local_var_count = frag.sg_blocks.len();
        let var_names = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j'];
        let renaming: HashMap<char, char> = (0..local_var_count)
            .map(|i| {
                (
                    var_names[i],
                    var_names[(self.var_offset + i) % var_names.len()],
                )
            })
            .collect();

        for sg in &frag.sg_blocks {
            self.sg_blocks
                .push(offset_sg_block(sg, self.atom_offset, &renaming));
        }
        for m in &frag.m_blocks {
            self.m_blocks.push(offset_m_block(m, self.atom_offset));
        }
        if let Some(c) = &frag.constraint {
            self.constraints.push(rename_constraint_vars(c, &renaming));
        }

        self.atom_offset += count_atoms(&frag.smiles);
        self.var_offset += local_var_count;
        let (text, rings_used) = offset_ring_labels(&frag.smiles, self.ring_offset);
        self.ring_offset += rings_used;
        self.smiles.push_str(&text);
    }

    /// `push_fragment` for a fragment that *is* one of the lipid's chains,
    /// additionally recording where its carbons landed. Fragments with no
    /// carbons (an empty slot's bare `O`) are appended but not recorded.
    fn push_chain(&mut self, frag: &CxFragment, sn: usize) {
        let base = self.atom_offset;
        let carbons: Vec<usize> = frag.carbon_atoms.iter().map(|a| a + base).collect();
        self.push_fragment(frag);
        if !carbons.is_empty() {
            // `base` is the slot's linking oxygen, or C1 for a chain pushed
            // without one. Either way it is the atom the chain hangs from,
            // which is what an sn position names.
            self.sn_sites.push((sn, base));
            self.chains.push((sn, carbons));
        }
    }

    /// Finalizes the assembly into the final string: base SMILES, the
    /// CXSMILES `|...|` block, then the trailing token list.
    ///
    /// The lipid-specific tokens sit in what a SMILES reader treats as the
    /// title field, and they only ever
    /// name things the `|...|` block also carries (`Sg:` variables,
    /// `$...$` atom labels), never atom positions, which no longer mean
    /// anything once a toolkit has renumbered the molecule.
    fn finish(self) -> Built {
        let mut blocks = Vec::new();
        let mut tokens = Vec::new();

        let mut sn_sites = self.sn_sites;
        sn_sites.sort();
        let swappable = self.swappable && sn_sites.len() >= 2;
        if swappable {
            blocks.push(sn_label_block(&sn_sites));
            tokens.push(format!(
                "swappable({})",
                sn_sites
                    .iter()
                    .map(|(sn, _)| format!("sn{sn}"))
                    .collect::<Vec<_>>()
                    .join(",")
            ));
        }

        blocks.extend(self.sg_blocks);
        blocks.extend(self.m_blocks);
        tokens.extend(self.constraints.iter().map(|c| format!("constrain({c})")));

        let mut chains = self.chains;
        chains.sort_by_key(|(sn, _)| *sn);

        let smiles = if blocks.is_empty() && tokens.is_empty() {
            self.smiles
        } else {
            format!(
                "{} |{}| {}",
                self.smiles,
                blocks.join(","),
                tokens.join(";")
            )
            .trim_end()
            .to_string()
        };
        Built {
            smiles,
            chains,
            sn_sites,
        }
    }
}

/// The `$...$` atom-label block naming each chain's linking atom `snN`.
///
/// Labels are positional and `;`-separated, one slot per atom of the main
/// SMILES in emission order; trailing empty slots are omitted, so the block
/// runs only as far as the last labelled atom. Unlike an atom *index*, a
/// label is maintained by the toolkits — canonicalizing the molecule
/// renumbers the atoms and rewrites this block to match, so `sn1` still
/// names the sn-1 position afterwards. That is why the `swappable` token
/// refers to these labels rather than to positions.
fn sn_label_block(sites: &[(usize, usize)]) -> String {
    let last = sites.iter().map(|(_, atom)| *atom).max().unwrap_or(0);
    let labels: Vec<String> = (0..=last)
        .map(|i| match sites.iter().find(|(_, atom)| *atom == i) {
            Some((sn, _)) => format!("sn{sn}"),
            None => String::new(),
        })
        .collect();
    format!("${}$", labels.join(";"))
}

/// Rewrites a `"Sg:n:LOCAL_POS:VAR:ht"` block to global atom numbering
/// and a possibly-renamed variable letter.
fn offset_sg_block(block: &str, atom_offset: usize, var_renaming: &HashMap<char, char>) -> String {
    let colons: Vec<usize> = block.match_indices(':').map(|(i, _)| i).collect();
    if colons.len() < 3 {
        return block.to_string();
    }
    let pos_start = colons[1] + 1;
    let pos_end = colons[2];
    let pos: usize = block[pos_start..pos_end].parse().unwrap_or(0);
    let var_start = pos_end + 1;
    let var_end = block[var_start..]
        .find(':')
        .map(|i| var_start + i)
        .unwrap_or(block.len());
    let var_char = block[var_start..var_end].chars().next().unwrap_or('a');
    let new_var = var_renaming.get(&var_char).copied().unwrap_or(var_char);
    format!("Sg:n:{}:{}:ht", pos + atom_offset, new_var)
}

/// Rewrites a `"m:LOCAL_ATOM_IDX:LOCAL_POS.LOCAL_POS..."` block to
/// global atom numbering.
fn offset_m_block(block: &str, atom_offset: usize) -> String {
    let colons: Vec<usize> = block.match_indices(':').map(|(i, _)| i).collect();
    if colons.len() < 2 {
        return block.to_string();
    }
    let idx_start = colons[0] + 1;
    let idx_end = colons[1];
    let atom_idx: usize = block[idx_start..idx_end].parse().unwrap_or(0);
    let positions_part = &block[idx_end + 1..];
    let adjusted_positions = positions_part
        .split('.')
        .map(|p| {
            p.parse::<usize>()
                .map(|v| (v + atom_offset).to_string())
                .unwrap_or_else(|_| p.to_string())
        })
        .collect::<Vec<_>>()
        .join(".");
    format!("m:{}:{}", atom_idx + atom_offset, adjusted_positions)
}

/// Shifts every `%nn` ring-closure label in a fragment up by `offset`,
/// returning the rewritten text and how many distinct rings it held.
///
/// A chain fragment always numbers its own rings from `%10` (see
/// `chain_tokens`), so a second ring-bearing chain in the same molecule
/// would otherwise reuse the first one's labels and bond the two chains
/// together. Rings are the only closures a chain emits, so shifting by the
/// running count is enough to keep them apart.
fn offset_ring_labels(smiles: &str, offset: usize) -> (String, usize) {
    if offset == 0 && !smiles.contains('%') {
        return (smiles.to_string(), 0);
    }
    let mut out = String::with_capacity(smiles.len());
    let mut seen: HashSet<usize> = HashSet::new();
    let mut rest = smiles;
    while let Some(at) = rest.find('%') {
        out.push_str(&rest[..at]);
        let digits = &rest[at + 1..];
        let end = digits
            .find(|c: char| !c.is_ascii_digit())
            .unwrap_or(digits.len());
        match digits[..end].parse::<usize>() {
            Ok(label) => {
                seen.insert(label);
                out.push_str(&format!("%{}", label + offset));
            }
            Err(_) => out.push('%'),
        }
        rest = &digits[end..];
    }
    out.push_str(rest);
    (out, seen.len())
}

/// Renames every variable letter in a `"a+b=N"`-style constraint per
/// `var_renaming`, leaving `+`/`=`/digits untouched.
fn rename_constraint_vars(constraint: &str, var_renaming: &HashMap<char, char>) -> String {
    constraint
        .chars()
        .map(|c| var_renaming.get(&c).copied().unwrap_or(c))
        .collect()
}

/// Every atom token in a SMILES string, across all dot-separated
/// components: each element symbol, each bracket atom, each `*` wildcard.
/// This is the indexing convention every `Sg:`/`m:`/`$...$` block uses.
///
/// Counting the leading character of each organic-subset symbol is enough
/// for what this crate emits: the two-letter halogens `Cl` and `Br` are
/// caught by their `C` and `B`, and neither trailing letter is itself in
/// the set. Ring-closure labels (`%10`) are punctuation and count for
/// nothing, which is what the `Sg:`/`m:` indices expect.
pub(crate) fn count_atoms(smiles: &str) -> usize {
    let mut count = 0;
    let mut chars = smiles.chars();

    while let Some(c) = chars.next() {
        if c == '[' {
            count += 1;
            for c2 in chars.by_ref() {
                if c2 == ']' {
                    break;
                }
            }
        } else if starts_atom(c) {
            count += 1;
        }
    }

    count
}

/// Whether `c` begins an atom token outside brackets.
///
/// The single definition of this convention; every caller that needs to
/// recognize atom-starting characters (including `insert_padding_atoms`)
/// goes through this function so the accepted character set can't drift
/// between call sites.
pub(crate) fn starts_atom(c: char) -> bool {
    matches!(
        c,
        'B' | 'C' | 'N' | 'O' | 'P' | 'S' | 'F' | 'I' | 'b' | 'c' | 'n' | 'o' | 's' | 'p' | '*'
    )
}

/// Resolves this generator's `Sg:n:pos:var:ht` unlocalized-double-bond
/// markers (see module docs and [`build_chain_segment`]) into one
/// concrete, fully-connected SMILES, for depiction purposes.
///
/// External SMILES renderers (chematic, RDKit, ...) have no notion of
/// these markers or of the `x+y=N` size constraint that goes with them —
/// parsing just the base SMILES and ignoring the CXSMILES suffix leaves
/// out most of the chain (the marker atoms stand in for a variable-length
/// run of `CH2` units whose *count* only exists in the constraint, not in
/// the base SMILES text), so the resulting depiction is a truncated stub.
///
/// Each Sg-marked atom is itself already the first unit of its own
/// variable's count (see `build_chain_segment`, which always emits
/// exactly one hardcoded atom per marker regardless of that variable's
/// eventual value), so only `value - 1` extra `C` atoms need inserting
/// after it. Any one valid split of each constraint sum across its
/// variables depicts an equally valid resolution, since the true
/// distribution is by definition unlocalized; this picks an even split.
///
/// An unlocalized modification's `m:` block is not resolved to a specific
/// carbon, so its dot-separated wildcard component remains unchanged.
///
/// Returns the input unchanged if it has no CXSMILES suffix or nothing
/// to expand.
pub fn smiles_expand(smi: &str) -> String {
    expand_with_padding_inserts(smi).0
}

/// [`smiles_expand`] plus the index translation it applied,
/// so callers holding atom indices into the unexpanded string can move them
/// across: the `(atom index, count)` runs of padding `C` atoms added for the
/// `Sg:` markers. These are chain carbons themselves, so a caller mapping a
/// chain also splices them into that chain's run (sort by position before
/// passing to [`shift_atom`]).
fn expand_with_padding_inserts(smi: &str) -> (String, Vec<(usize, usize)>) {
    let Some(suffix_start) = smi.find(" |") else {
        return (smi.to_string(), Vec::new());
    };
    let Some(rest) = smi.get(suffix_start + 2..) else {
        return (smi[..suffix_start].to_string(), Vec::new());
    };
    let Some(close) = rest.find('|') else {
        return (smi[..suffix_start].to_string(), Vec::new());
    };
    let base = &smi[..suffix_start];
    let blocks = &rest[..close];
    let constraints = &rest[close + 1..];

    let sg_positions: Vec<usize> = blocks
        .split(',')
        .filter_map(|b| b.strip_prefix("Sg:n:"))
        .filter_map(|r| r.split(':').next())
        .filter_map(|pos| pos.parse::<usize>().ok())
        .collect();

    if sg_positions.is_empty() {
        return (base.to_string(), Vec::new());
    }

    // Each `constrain(a+b=N)` corresponds, in order, to the next `n_terms`
    // Sg positions: variable *names* can repeat across independently
    // numbered chains (`CxBuilder`'s var_offset wraps at ten), so positional
    // matching against emission order — not name matching — is what stays
    // unambiguous here.
    let mut inserts: Vec<(usize, usize)> = Vec::new();
    let mut consumed = 0usize;
    for (n_terms, total) in trailer_equations(constraints) {
        if n_terms == 0 || consumed + n_terms > sg_positions.len() {
            break;
        }
        let each = total / n_terms;
        let rem = total % n_terms;
        for k in 0..n_terms {
            let value = each + if k + 1 == n_terms { rem } else { 0 };
            inserts.push((sg_positions[consumed + k], value.saturating_sub(1)));
        }
        consumed += n_terms;
    }

    let padded = if inserts.is_empty() {
        base.to_string()
    } else {
        insert_padding_atoms(base, &inserts)
    };
    (padded, inserts)
}

/// The `constrain(...)` equations of a trailer, as `(term count, sum)` in
/// order.
///
/// The trailer is normally a `;`-separated token list. Comma-separated bare
/// `a+b=15` equations are accepted as a compatibility spelling.
pub(crate) fn trailer_equations(trailer: &str) -> Vec<(usize, usize)> {
    let equation = |eq: &str| -> Option<(usize, usize)> {
        let (lhs, rhs) = eq.split_once('=')?;
        Some((lhs.split('+').count(), rhs.trim().parse().ok()?))
    };
    let mut out = Vec::new();
    for token in trailer.trim().split(';') {
        let token = token.trim();
        if token.is_empty() {
            continue;
        }
        match token.split_once('(') {
            Some(("constrain", body)) => out.extend(equation(body.trim_end_matches(')'))),
            Some(_) => {} // some other token; unknown tokens are ignored by design
            None => out.extend(token.split(',').filter_map(equation)),
        }
    }
    out
}

/// Whether a trailer carries a `swappable(...)` token, i.e. whether the sn
/// assignment the string shows is one arbitrary choice among several.
pub(crate) fn trailer_is_swappable(trailer: &str) -> bool {
    trailer
        .split(';')
        .any(|t| t.trim().starts_with("swappable("))
}

/// Translates a pre-expansion atom index. `shift_table` is the `(position,
/// count)` padding runs sorted by position; `insert_padding_atoms` writes
/// each run *after* the atom at `position`, so an index only shifts by the
/// padding that precedes it.
fn shift_atom(atom: usize, shift_table: &[(usize, usize)]) -> usize {
    atom + shift_table
        .iter()
        .take_while(|(pos, _)| *pos < atom)
        .map(|(_, count)| *count)
        .sum::<usize>()
}

/// Inserts `count` extra `C` atoms immediately after atom index `pos`
/// (0-based, counting every element/bracket-atom/`*` token across the
/// whole SMILES — same convention as `count_atoms`) for
/// each `(pos, count)` pair.
fn insert_padding_atoms(smi: &str, inserts: &[(usize, usize)]) -> String {
    let insert_map: HashMap<usize, usize> = inserts.iter().copied().collect();
    let chars: Vec<char> = smi.chars().collect();
    let padding_total: usize = inserts.iter().map(|(_, n)| n).sum();
    let mut out = String::with_capacity(smi.len() + padding_total);
    let mut atom_idx = 0usize;
    let mut i = 0usize;
    while i < chars.len() {
        let c = chars[i];
        if c == '[' {
            out.push(c);
            i += 1;
            while i < chars.len() {
                out.push(chars[i]);
                let closed = chars[i] == ']';
                i += 1;
                if closed {
                    break;
                }
            }
            if let Some(&n) = insert_map.get(&atom_idx) {
                out.extend(std::iter::repeat_n('C', n));
            }
            atom_idx += 1;
            continue;
        }
        if starts_atom(c) {
            out.push(c);
            if let Some(&n) = insert_map.get(&atom_idx) {
                out.extend(std::iter::repeat_n('C', n));
            }
            atom_idx += 1;
        } else {
            out.push(c);
        }
        i += 1;
    }
    out
}

// ---------- chain-token parsing ----------

fn parse_position_list(s: &str) -> Option<Vec<(u32, Option<Geometry>)>> {
    if s.is_empty() {
        return Some(Vec::new());
    }
    s.split(',')
        .map(|tok| {
            let tok = tok.trim();
            let (num_part, geom) = if let Some(stripped) = tok.strip_suffix(['Z', 'z']) {
                (stripped, Some(Geometry::Cis))
            } else if let Some(stripped) = tok.strip_suffix(['E', 'e']) {
                (stripped, Some(Geometry::Trans))
            } else {
                (tok, None)
            };
            num_part.parse::<u32>().ok().map(|n| (n, geom))
        })
        .collect()
}

fn parse_plain_position_list(s: &str) -> Option<Vec<u32>> {
    if s.is_empty() {
        return Some(Vec::new());
    }
    s.split(',').map(|t| t.trim().parse::<u32>().ok()).collect()
}

/// Fill in `missing` additional double-bond positions beyond whatever is
/// already in `existing`, avoiding any existing position or its immediate
/// neighbors (no cumulated dienes). Starts at C9 and steps by 3 (the
/// common methylene-interrupted spacing), wrapping back to the start of
/// the valid range if it runs off the end of the chain. Returns fewer
/// than `missing` positions only if the chain is too short/crowded to fit
/// them all.
fn placeholder_positions(carbon: u32, existing: &[u32], missing: usize) -> Vec<u32> {
    if missing == 0 || carbon < 4 {
        return Vec::new();
    }
    let max_pos = carbon.saturating_sub(1);
    let mut used: HashSet<u32> = existing.iter().copied().collect();
    let mut result = Vec::new();
    let mut next_natural: u32 = 9.min(max_pos).max(2);

    while result.len() < missing {
        let mut candidate = next_natural;
        let mut found = None;
        for _ in 0..=(2 * max_pos.max(2)) {
            if candidate < 2 || candidate > max_pos {
                candidate = 2;
            }
            let free = !used.contains(&candidate)
                && !used.contains(&candidate.saturating_sub(1))
                && !used.contains(&(candidate + 1));
            if free {
                found = Some(candidate);
                break;
            }
            candidate += 1;
        }
        match found {
            Some(p) => {
                used.insert(p);
                result.push(p);
                next_natural = p + 3;
            }
            None => break,
        }
    }
    result
}

fn parse_chain_token(tok: &str) -> Option<ParsedChain> {
    let tok = tok.trim();
    if tok == "0:0" {
        return Some(ParsedChain {
            prefix: ChainPrefix::Acyl,
            carbon: 0,
            db_pos: Vec::new(),
            mods: Vec::new(),
            rings: Vec::new(),
        });
    }

    let mut prefix = ChainPrefix::Acyl;
    let mut s = tok;
    if let Some(rest) = s.strip_prefix("O-") {
        prefix = ChainPrefix::EtherAlkyl;
        s = rest;
    } else if let Some(rest) = s.strip_prefix("P-") {
        prefix = ChainPrefix::EtherAlkenyl;
        s = rest;
    } else if (s.starts_with('d') || s.starts_with('D'))
        && s[1..].chars().next().is_some_and(|c| c.is_ascii_digit())
    {
        prefix = ChainPrefix::SphingoidD;
        s = &s[1..];
    } else if (s.starts_with('t') || s.starts_with('T'))
        && s[1..].chars().next().is_some_and(|c| c.is_ascii_digit())
    {
        prefix = ChainPrefix::SphingoidT;
        s = &s[1..];
    }

    let (c_str, rest) = s.split_once(':')?;
    let carbon: u32 = c_str.parse().ok()?;

    let db_digits_end = rest
        .find(|c: char| !c.is_ascii_digit())
        .unwrap_or(rest.len());
    let db_declared: u32 = rest[..db_digits_end].parse().ok()?;
    let mut cursor = &rest[db_digits_end..];

    let mut given_db = Vec::new();
    if let Some(after_paren) = cursor.strip_prefix('(') {
        let end = after_paren.find(')')?;
        given_db = parse_position_list(&after_paren[..end])?;
        cursor = &after_paren[end + 1..];
    }

    let mut mods = Vec::new();
    let mut rings = Vec::new();
    let mut ring_db = Vec::new();
    let mut oxy_declared: u32 = 0;

    for seg in cursor.split(';') {
        // A ring's structural detail is written in square brackets with the
        // groups sitting on it inside — `FA 20:2;[8-12cy5;11OH;9oxo];15OH`.
        // Those inner groups are ordinary position/abbreviation pairs, so
        // dropping the brackets after the `;` split is enough to read them.
        let seg = seg.trim().trim_start_matches('[').trim_end_matches(']');
        if seg.is_empty() {
            continue;
        }
        parse_modification_segment(
            seg,
            carbon,
            &mut mods,
            &mut rings,
            &mut ring_db,
            &mut oxy_declared,
        )?;
    }

    if given_db.len() > db_declared as usize {
        return None; // more positions given than the declared count -- malformed
    }
    let missing = db_declared as usize - given_db.len();
    let existing_positions: Vec<u32> = given_db.iter().map(|&(p, _)| p).collect();
    let extra = placeholder_positions(carbon, &existing_positions, missing);
    if extra.len() < missing {
        return None; // chain too short/crowded to even place a placeholder
    }
    let mut db_pos: Vec<DbPos> = given_db
        .into_iter()
        .map(|(pos, geom)| DbPos {
            pos,
            geom,
            placeholder: false,
        })
        .collect();
    db_pos.extend(extra.into_iter().map(|pos| DbPos {
        pos,
        geom: None,
        placeholder: true,
    }));
    // A ring's own double bonds sit outside the chain's `C:DBE` count, so
    // they join only after that count has been reconciled.
    db_pos.extend(ring_db);

    // A generic `;O2` says only how many oxygens there are. That is a
    // hypothesis-free statement unless some group was also named, whether
    // positioned or not.
    if oxy_declared > 0 && mods.is_empty() && rings.is_empty() {
        return None; // generic oxygen count with no position hypothesis
    }

    Some(ParsedChain {
        prefix,
        carbon,
        db_pos,
        mods,
        rings,
    })
}

/// Reads one `;`-separated functional-group segment into `mods`/`rings`.
///
/// Recognized forms, in the order they are tried:
///
/// * `8-12cy5`, `11-13cy3:0`, `9-11cy3:1(9)` — a Table 1B carbocycle
///   spanning those chain carbons, optionally with double bonds of its own.
/// * `(OH)2`, `(NO2)` — the paper's parenthesized form, used when a group
///   occurs more than once (followed by the count) or when its abbreviation
///   itself contains digits. Positions are undetermined.
/// * `OH(3,5)`, `ep(5)` — compatibility aliases for positioned groups.
/// * `11OH,15OH`, `3Me,7Me`, `oxo` — position (optional) in front of the
///   abbreviation, comma-separated, which is what the paper specifies.
/// * `O`, `O2` — a bare count of oxygens with no group breakdown at all.
///
/// Returns `None` for a malformed segment, and for the Table 1A groups
/// listed in [`UNRENDERABLE`].
fn parse_modification_segment(
    seg: &str,
    carbon: u32,
    mods: &mut Vec<Mod>,
    rings: &mut Vec<Ring>,
    ring_db: &mut Vec<DbPos>,
    oxy_declared: &mut u32,
) -> Option<()> {
    if let Some(ring) = parse_ring_segment(seg, carbon, ring_db)? {
        rings.push(ring);
        return Some(());
    }

    if let Some(rest) = seg.strip_prefix('(') {
        let (abbr, count) = rest.split_once(')')?;
        let count: u32 = if count.is_empty() {
            1
        } else {
            count.parse().ok()?
        };
        let branch = substituent(abbr)?;
        mods.extend((0..count).map(|_| Mod { pos: 0, branch }));
        return Some(());
    }

    if let Some((abbr, rest)) = seg.split_once('(') {
        let positions = parse_plain_position_list(rest.strip_suffix(')')?)?;
        if let Some(ring_size) = alias_ring_size(abbr) {
            for p in positions {
                rings.push(ring_at(p, ring_size, carbon)?);
            }
        } else {
            let branch = substituent(abbr)?;
            mods.extend(positions.into_iter().map(|pos| Mod { pos, branch }));
        }
        return Some(());
    }

    for part in seg.split(',') {
        let digits_end = part.find(|c: char| !c.is_ascii_digit()).unwrap_or(0);
        let (pos, abbr) = part.split_at(digits_end);
        let pos: u32 = if pos.is_empty() { 0 } else { pos.parse().ok()? };
        if let Some(ring_size) = alias_ring_size(abbr) {
            rings.push(ring_at(pos, ring_size, carbon)?);
        } else if let Some(branch) = substituent(abbr) {
            mods.push(Mod { pos, branch });
        } else if UNRENDERABLE.contains(&abbr) {
            return None;
        } else if let Some(digits) = abbr.strip_prefix('O') {
            // A bare oxygen count, e.g. `;O2`. It carries no position, and
            // `pos` is part of the count's own token, so nothing is stored.
            *oxy_declared = if digits.is_empty() {
                1
            } else {
                digits.parse().ok()?
            };
        } else {
            return None; // unrecognized functional group
        }
    }
    Some(())
}

/// The `SUBSTITUENTS` branch for a Table 1A abbreviation or accepted alias.
fn substituent(abbr: &str) -> Option<&'static str> {
    match abbr {
        "hydroxy" => Some("O"),
        "keto" => Some("=O"),
        _ => SUBSTITUENTS
            .iter()
            .find(|(name, _)| *name == abbr)
            .map(|(_, branch)| *branch),
    }
}

/// Ring size for the `ep(5)`/`cyc(5)` compatibility aliases.
fn alias_ring_size(abbr: &str) -> Option<u32> {
    match abbr {
        "Ep" | "ep" | "epox" => Some(2),
        "cyc" | "cyclo" => Some(3),
        _ => None,
    }
}

/// A ring of `size` carbons beginning at `start`. Size 2 means an epoxide:
/// two carbons bridged by an oxygen.
fn ring_at(start: u32, size: u32, carbon: u32) -> Option<Ring> {
    let end = start + size - 1;
    if start == 0 || end > carbon {
        return None; // ring runs off the end of the chain
    }
    Some(Ring {
        start,
        end,
        bridge: (size == 2).then_some("O"),
    })
}

/// A Table 1B `cyX` segment: `<start>-<end>cy<ring atoms>[:<DBEs>[(<pos>)]]`,
/// e.g. `11-13cy3:0` (lactobacillic acid) or `9-11cy3:1(9)` (sterculic).
/// Returns `Some(None)` when the segment is not a ring at all.
///
/// A ring's own double bonds are declared inside it and are *additional* to
/// the chain's `C:DBE` count, so they are collected separately in `ring_db`
/// and merged after that count has been checked.
fn parse_ring_segment(seg: &str, carbon: u32, ring_db: &mut Vec<DbPos>) -> Option<Option<Ring>> {
    let Some((span, rest)) = seg.split_once("cy") else {
        return Some(None);
    };
    let Some((start, end)) = span.split_once('-') else {
        return Some(None);
    };
    let start: u32 = start.trim().parse().ok()?;
    let end: u32 = end.trim().parse().ok()?;

    let (size, dbe) = match rest.split_once(':') {
        Some((size, dbe)) => (size, dbe),
        None => (rest, ""),
    };
    // `cy5` states the ring's atom count; the span has to agree with it.
    if end < start || end > carbon || end - start + 1 != size.trim().parse::<u32>().ok()? {
        return None;
    }

    let (declared, positions) = match dbe.split_once('(') {
        Some((n, pos)) => (n, parse_position_list(pos.strip_suffix(')')?)?),
        None => (dbe, Vec::new()),
    };
    let declared: usize = if declared.trim().is_empty() {
        0
    } else {
        declared.trim().parse().ok()?
    };
    if positions.len() > declared {
        return None;
    }
    // An unpositioned ring double bond has only one place it can sit that
    // the name does not already pin down, so it is not a guess to leave the
    // count unmatched — but it is also not localizable, so it is dropped
    // rather than invented. Only stated positions are rendered.
    for (pos, geom) in positions {
        if pos < start || pos >= end {
            return None; // a ring's double bond has to be inside the ring
        }
        ring_db.push(DbPos {
            pos,
            geom,
            placeholder: false,
        });
    }
    Some(Some(Ring {
        start,
        end,
        bridge: None,
    }))
}

// ---------- generic chain-body rendering ----------

/// Per-carbon SMILES atom tokens (1-indexed via `atoms[k-1]`) plus the
/// bond text keyed by the *starting* carbon of each bond (`k` for the
/// bond between Ck and Ck+1); missing keys mean a plain single bond.
/// `db_pos` must only contain known (non-placeholder) double bonds —
/// unlocalized ones are represented separately via `Sg:` blocks (see
/// `build_chain_segment`), never as a literal position here.
///
/// A carbon carrying functional groups gets them as branches in table
/// order, after any ring-closure label: `C%10(O)(C)`. Ring labels use the
/// two-digit `%nn` form starting at `%10` because the single digits are
/// already spoken for by the fixed headgroup templates (the sterol nucleus
/// alone uses 1, 2 and 3); `CxBuilder::push_fragment` renumbers them
/// upward as fragments are assembled, exactly as it does the `Sg:`
/// variable letters.
fn chain_tokens(
    carbon: u32,
    db_pos: &[DbPos],
    mods: &[Mod],
    rings: &[Ring],
) -> (Vec<String>, HashMap<u32, String>) {
    let db_start: HashMap<u32, DbPos> = db_pos.iter().map(|d| (d.pos, *d)).collect();

    let mut labels: HashMap<u32, String> = HashMap::new();
    let mut bonds: HashMap<u32, String> = HashMap::new();
    for (i, ring) in rings.iter().enumerate() {
        let label = format!("%{}", 10 + i);
        labels.entry(ring.start).or_default().push_str(&label);
        labels.entry(ring.end).or_default().push_str(&label);
        // An epoxide's two carbons are adjacent, so its oxygen goes in the
        // chain bond between them and the ring closure supplies the C-C
        // bond: `C%10OC%10`. A carbocycle needs no bridge — the chain
        // itself is the rest of the ring.
        if let Some(bridge) = ring.bridge {
            bonds.insert(ring.start, bridge.to_string());
        }
    }

    let mut atoms = Vec::with_capacity(carbon as usize);
    for k in 1..=carbon {
        let mut tok = String::from("C");
        if let Some(label) = labels.get(&k) {
            tok.push_str(label);
        }
        for m in mods.iter().filter(|m| m.pos == k) {
            tok.push('(');
            tok.push_str(m.branch);
            tok.push(')');
        }
        atoms.push(tok);
    }

    for &k in db_start.keys() {
        bonds.insert(k, "=".to_string());
    }
    for (&k, d) in &db_start {
        if let Some(g) = d.geom {
            if k > 1 {
                bonds.entry(k - 1).or_insert_with(|| "/".to_string());
            }
            if k + 1 < carbon {
                let sym = match g {
                    Geometry::Cis => "\\",
                    Geometry::Trans => "/",
                };
                bonds.entry(k + 1).or_insert_with(|| sym.to_string());
            }
        }
    }
    (atoms, bonds)
}

/// Emits `start..=end` and, alongside, the local atom index each carbon
/// lands on. A carbon's own token always leads with its `C`, so the index
/// is simply the atom count of everything emitted before it; any branch
/// atoms the token carries (`(=O)`, `(O)`, ...) follow and are skipped over
/// by the running count, as does an epoxide's bridging oxygen when the
/// bond text between two carbons carries an atom of its own.
fn assemble_range(
    atoms: &[String],
    bonds: &HashMap<u32, String>,
    start: u32,
    end: u32,
) -> (String, Vec<usize>) {
    let mut out = String::new();
    let mut carbon_atoms = Vec::with_capacity((end.saturating_sub(start) + 1) as usize);
    let mut atom_idx = 0usize;
    for k in start..=end {
        let tok = &atoms[(k - 1) as usize];
        carbon_atoms.push(atom_idx);
        atom_idx += count_atoms(tok);
        out.push_str(tok);
        if k < end {
            if let Some(b) = bonds.get(&k) {
                out.push_str(b);
                atom_idx += count_atoms(b);
            }
        }
    }
    (out, carbon_atoms)
}

/// Builds the SMILES for the `start..=end` slice of a chain's `1..=carbon`
/// numbering (bond/geometry computation still spans the whole chain, so a
/// `start` past 1 stays bond-consistent with whatever precedes it).
/// `carbonyl_c1` makes C1 a carbonyl carbon (acyl linkage via
/// ester/amide) when `start == 1`; otherwise C1 is a plain sp3 carbon
/// (ether linkage) or, for `start > 1`, irrelevant (C1 isn't emitted).
/// Used directly with `start == 1, end == carbon` for a normal chain, and
/// with `start == 3` for a sphingoid base's own tail (C1/C2 are handled
/// by the caller's fixed template).
#[allow(clippy::too_many_arguments)]
fn build_chain_range(
    carbon: u32,
    db_pos: &[DbPos],
    mods: &[Mod],
    rings: &[Ring],
    carbonyl_c1: bool,
    start: u32,
    end: u32,
) -> (String, Vec<usize>) {
    if carbon == 0 || end < start {
        return (String::new(), Vec::new());
    }
    let (mut atoms, bonds) = chain_tokens(carbon, db_pos, mods, rings);
    if carbonyl_c1 {
        atoms[0] = "C(=O)".to_string();
    }
    assemble_range(&atoms, &bonds, start, end)
}

// ---------- headgroup builders ----------
//
// Each builder writes one lipid class's fixed template and delegates its
// variable slots to `build_chain_slot`/`build_chain_fragment`, which return a
// chain fragment starting at C1 ready to be attached after the linking
// heteroatom (`O` for ester/ether, `N` for amide) the template itself writes.

fn fa_smiles(chains: &[ParsedChain]) -> Option<Built> {
    let c = chains.first()?;
    if c.carbon == 0 || c.prefix != ChainPrefix::Acyl {
        return None;
    }
    let frag = build_chain_fragment(c)?;
    let mut b = CxBuilder::default();
    b.push_fixed("O");
    b.push_chain(&frag, 1);
    Some(b.finish())
}

/// AMPP (N-(4-aminomethylphenyl)pyridinium) charge-tagged fatty acid, the
/// derivatization used by the `hete-ead`/`ara-ead` EAD test fixtures
/// (`AMP-FA 16:0` in `data/AMP_PUFA.csv`). The acyl chain forms an amide
/// with the tag's benzylamine; the pyridinium ring carries the fixed
/// positive charge.
fn amp_fa_smiles(chains: &[ParsedChain]) -> Option<Built> {
    let c = chains.first()?;
    if c.carbon == 0 || c.prefix != ChainPrefix::Acyl {
        return None;
    }
    let frag = build_chain_fragment(c)?;
    let mut b = CxBuilder::default();
    b.push_fixed("[n+]1ccccc1-c1ccc(CN");
    b.push_chain(&frag, 1);
    b.push_fixed(")cc1");
    Some(b.finish())
}

fn nae_smiles(chains: &[ParsedChain]) -> Option<Built> {
    let c = chains.first()?;
    if c.carbon == 0 || c.prefix != ChainPrefix::Acyl {
        return None;
    }
    let frag = build_chain_fragment(c)?;
    let mut b = CxBuilder::default();
    b.push_fixed("OCCN");
    b.push_chain(&frag, 1);
    Some(b.finish())
}

fn car_smiles(chains: &[ParsedChain]) -> Option<Built> {
    let c = chains.first()?;
    if c.carbon == 0 || c.prefix != ChainPrefix::Acyl {
        return None;
    }
    let frag = build_chain_fragment(c)?;
    // The ester O needs two real bonds: to the chain's carbonyl carbon
    // (C1, i.e. frag's own first atom) and to carnitine's chiral carbon.
    // Writing the chain as a branch straight off O keeps both bonds
    // explicit without touching [C@H]'s neighbor order (O still precedes
    // it in the text exactly as before), so no stereo re-derivation is
    // needed.
    let mut b = CxBuilder::default();
    b.push_fixed("O(");
    b.push_chain(&frag, 1);
    b.push_fixed(")[C@H](CC(=O)[O-])C[N+](C)(C)C");
    Some(b.finish())
}

fn ce_smiles(chains: &[ParsedChain]) -> Option<Built> {
    let c = chains.first()?;
    if c.carbon == 0 || c.prefix != ChainPrefix::Acyl {
        return None;
    }
    let frag = build_chain_fragment(c)?;
    let mut b = CxBuilder::default();
    b.push_fixed("C12(CC=C3CC(O");
    b.push_chain(&frag, 1);
    b.push_fixed(")CCC3(C)C1CCC1(C)C(C(C)CCCC(C)C)CCC21)");
    Some(b.finish())
}

fn st_smiles() -> String {
    "C12(CC=C3CC(O)CCC3(C)C1CCC1(C)C(C(C)CCCC(C)C)CCC21)".to_string()
}

/// Glycerol backbone with up to `slots` ester/ether positions (sn1, sn2,
/// sn3 in that order); unfilled positions are free hydroxyls. Used for
/// MG/DG/TG, which have no phosphate headgroup.
///
/// `swappable` records that the name joined its chains with `_`, so the sn
/// order shown is one arbitrary choice. There is only one layout either
/// way: the chains are written into the main string, where their `Sg:`/`m:`
/// blocks are legal, and the ambiguity is stated by labelling the linking
/// atoms and emitting a `swappable(...)` token.
fn glycerolipid_smiles(chains: &[ParsedChain], slots: usize, swappable: bool) -> Option<Built> {
    if chains.is_empty() || chains.len() > slots {
        return None;
    }
    // A slot beyond this class's count is a free hydroxyl, which is
    // exactly what `build_chain_slot(None)` builds.
    let sn1 = build_chain_slot(chains.first())?;
    let sn2 = build_chain_slot(chains.get(1).filter(|_| slots >= 2))?;
    let sn3 = build_chain_slot(chains.get(2).filter(|_| slots >= 3))?;

    let mut b = CxBuilder {
        swappable,
        ..Default::default()
    };
    b.push_fixed("C(C");
    b.push_chain(&sn3, 3);
    b.push_fixed(")(");
    b.push_chain(&sn2, 2);
    b.push_fixed(")C");
    b.push_chain(&sn1, 1);
    Some(b.finish())
}

/// Diacyl-glycerophospholipid: `sn1`/`sn2` chains plus a fixed
/// phospho-headgroup tail attached at sn3 (e.g. `"OP(=O)([O-])OCC[N+](C)(C)C"`
/// for PC). Also covers the lyso forms, which supply one chain and so have
/// nothing to be ambiguous between. `headgroup_tail` is generic across
/// PC/PE/PS/PG/PI/PA and their lyso forms, so this one builder gives all of
/// them unlocalized-double-bond and regiochemistry coverage.
fn gpl_smiles(headgroup_tail: &str, chains: &[ParsedChain], swappable: bool) -> Option<Built> {
    if chains.is_empty() || chains.len() > 2 {
        return None;
    }
    let sn1 = build_chain_slot(chains.first())?;
    let sn2 = build_chain_slot(chains.get(1))?;
    let mut b = CxBuilder {
        swappable,
        ..Default::default()
    };
    b.push_fixed("C(C");
    b.push_fixed(headgroup_tail);
    b.push_fixed(")(");
    b.push_chain(&sn2, 2);
    b.push_fixed(")C");
    b.push_chain(&sn1, 1);
    Some(b.finish())
}

/// Cardiolipin: two phosphatidyl arms (sn1/sn2 and sn3/sn4) hung off a
/// central glycerol's C1/C3, with a free hydroxyl at the central C2.
fn cl_smiles(chains: &[ParsedChain], swappable: bool) -> Option<Built> {
    if chains.len() != 4 {
        return None;
    }
    let sn1 = build_chain_slot(chains.first())?;
    let sn2 = build_chain_slot(chains.get(1))?;
    let sn3 = build_chain_slot(chains.get(2))?;
    let sn4 = build_chain_slot(chains.get(3))?;

    let mut b = CxBuilder {
        swappable,
        ..Default::default()
    };
    b.push_fixed("C(COP(=O)(O)OCC(");
    b.push_chain(&sn2, 2);
    b.push_fixed(")C");
    b.push_chain(&sn1, 1);
    b.push_fixed(")(O)COP(=O)(O)OCC(");
    b.push_chain(&sn4, 4);
    b.push_fixed(")C");
    b.push_chain(&sn3, 3);
    Some(b.finish())
}

/// Sphingoid-base backbone (C1..Cn), with the N-acyl amide (or, for free
/// bases, a plain amine) at C2 and a fixed hydroxyl at C3 (plus C4 for
/// trihydroxy `t`-prefixed bases), plus whatever additional db/OH/oxo the
/// base chain itself declares beyond C3/C4. `head_tail` is appended after
/// C1's own linking oxygen (e.g. `"O"` for a free primary alcohol,
/// `"OP(=O)([O-])OCC[N+](C)(C)C"` for sphingomyelin's phosphocholine).
/// The sphingoid base and its N-acyl chain are chemically distinct
/// positions (never ambiguous with each other), so this has no
/// regiochemistry mode. Both the base's own tail (C3 onward) and the
/// N-acyl chain are Sg:-aware.
fn sphingoid_smiles(
    base: &ParsedChain,
    n_acyl: Option<&ParsedChain>,
    head_tail: &str,
) -> Option<Built> {
    if base.carbon < 4 {
        return None;
    }
    let is_triol = base.prefix == ChainPrefix::SphingoidT;
    let mut mods = base.mods.clone();
    for fixed in [3, 4] {
        // C3 always carries a hydroxyl; C4 does too on a `t`-prefixed triol.
        if (fixed == 3 || is_triol) && !mods.iter().any(|m| m.pos == fixed) {
            mods.push(Mod {
                pos: fixed,
                branch: "O",
            });
        }
    }

    let rest = build_chain_segment(3, base.carbon, &base.db_pos, &mods, &base.rings, false)?;

    let n_frag = match n_acyl {
        Some(acyl) if acyl.carbon > 0 && acyl.prefix == ChainPrefix::Acyl => {
            Some(build_chain_fragment(acyl)?)
        }
        Some(_) => return None,
        None => None,
    };

    let mut b = CxBuilder::default();
    // The template writes C2 first (it carries the N branch), then C1; the
    // base's own `rest` fragment picks up from C3, so the base chain's
    // C1..Cn run starts with those two fixed atoms.
    b.push_fixed("C(C");
    b.push_fixed(head_tail);
    b.push_fixed(")(N");
    if let Some(f) = &n_frag {
        b.push_chain(f, 2);
    }
    b.push_fixed(")");
    let rest_offset = b.atom_offset;
    b.push_fragment(&rest);
    let mut base_carbons = vec![1usize, 0usize];
    base_carbons.extend(rest.carbon_atoms.iter().map(|a| a + rest_offset));
    b.chains.push((1, base_carbons));
    Some(b.finish())
}

fn sphingo_dispatch(chains: &[ParsedChain], head_tail: &str, has_n_acyl: bool) -> Option<Built> {
    let base = chains.first()?;
    if base.carbon == 0 {
        return None;
    }
    let n_acyl = if has_n_acyl { chains.get(1) } else { None };
    if has_n_acyl && n_acyl.is_none() {
        return None;
    }
    sphingoid_smiles(base, n_acyl, head_tail)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Cheap structural sanity check (not a real SMILES parser): every
    /// bracket type balances, and every numeric ring-closure digit is
    /// used an even number of times.
    fn assert_balanced(smiles: &str) {
        let body = smiles.split(" |").next().unwrap_or(smiles);
        let mut paren = 0i32;
        let mut bracket = 0i32;
        for c in body.chars() {
            match c {
                '(' => paren += 1,
                ')' => paren -= 1,
                '[' => bracket += 1,
                ']' => bracket -= 1,
                _ => {}
            }
            assert!(paren >= 0, "unbalanced ) in {smiles}");
            assert!(bracket >= 0, "unbalanced ] in {smiles}");
        }
        assert_eq!(paren, 0, "unbalanced parens in {smiles}");
        assert_eq!(bracket, 0, "unbalanced brackets in {smiles}");

        let mut ring_counts: HashMap<char, u32> = HashMap::new();
        for c in body.chars().filter(|c| c.is_ascii_digit()) {
            *ring_counts.entry(c).or_insert(0) += 1;
        }
        for (digit, n) in ring_counts {
            assert_eq!(
                n % 2,
                0,
                "ring closure digit {digit} used {n} times in {smiles}"
            );
        }
    }

    /// The element symbol of every atom in a plain SMILES, in emission
    /// order — the indexing `LipidStructure::chains` refers to.
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
                    out.push(
                        chars[start..=i.min(chars.len() - 1)]
                            .iter()
                            .collect::<String>(),
                    );
                    i += 1;
                }
                c @ ('B' | 'C' | 'N' | 'O' | 'P' | 'S' | 'F' | 'I' | 'c' | 'n' | 'o' | 's'
                | 'p' | '*') => {
                    out.push(c.to_string());
                    i += 1;
                }
                _ => i += 1,
            }
        }
        out
    }

    /// Every mapped carbon must actually be a carbon, the run must cover
    /// the whole chain exactly once, and chains must not overlap.
    fn assert_chain_map_sane(name: &str, expected: &[(usize, usize)]) {
        let st = generate_structure(name).unwrap_or_else(|| panic!("{name} should resolve"));
        assert_balanced(&st.smiles);
        assert!(
            !st.smiles.contains(" |"),
            "{name} left a CXSMILES suffix: {}",
            st.smiles
        );
        let symbols = atom_symbols(&st.smiles);
        assert_eq!(
            symbols.len(),
            st.atom_count,
            "{name} atom_count disagrees with the SMILES"
        );

        assert_eq!(st.chains.len(), expected.len(), "{name} chain count");
        let mut seen: HashSet<usize> = HashSet::new();
        for (chain, &(sn, carbons)) in st.chains.iter().zip(expected) {
            assert_eq!(chain.sn, sn, "{name} sn position");
            assert_eq!(chain.carbons.len(), carbons, "{name} sn{sn} carbon count");
            for (k, &atom) in chain.carbons.iter().enumerate() {
                assert!(
                    atom < symbols.len(),
                    "{name} sn{sn} C{} out of range",
                    k + 1
                );
                assert_eq!(
                    symbols[atom],
                    "C",
                    "{name} sn{sn} C{} maps to atom {atom} ({}), not a carbon, in {}",
                    k + 1,
                    symbols[atom],
                    st.smiles
                );
                assert!(
                    seen.insert(atom),
                    "{name} atom {atom} claimed by two chains"
                );
            }
        }
    }

    #[test]
    fn structure_maps_single_acyl_chain() {
        assert_chain_map_sane("FA 18:1(9)", &[(1, 18)]);
    }

    #[test]
    fn structure_maps_both_gpl_chains_in_sn_order() {
        assert_chain_map_sane("PC 16:0/18:1(9)", &[(1, 16), (2, 18)]);
    }

    #[test]
    fn structure_maps_all_three_glycerolipid_chains() {
        assert_chain_map_sane("TG 16:0/18:1(9)/18:2(9,12)", &[(1, 16), (2, 18), (3, 18)]);
    }

    #[test]
    fn structure_maps_sphingoid_base_and_n_acyl() {
        assert_chain_map_sane("Cer d18:1(4)/16:0", &[(1, 18), (2, 16)]);
    }

    #[test]
    fn structure_maps_chain_with_modifications() {
        assert_chain_map_sane("AMP-FA 20:4(5,8,11,14);15OH", &[(1, 20)]);
    }

    /// An unlocalized chain is emitted as a short `Sg:`-marked stub; the
    /// map has to pick up the padding atoms that expansion splices in, or
    /// the chain would come back several carbons short.
    #[test]
    fn structure_map_covers_padding_atoms_from_expansion() {
        assert_chain_map_sane("FA 18:1", &[(1, 18)]);
        assert_chain_map_sane("PC 16:0/18:2", &[(1, 16), (2, 18)]);
    }

    /// The structure view drops the sn labels and the token, keeping only
    /// the connected molecule, and says the assignment is arbitrary via
    /// `regio_resolved` instead.
    #[test]
    fn structure_builds_sn_unresolved_names_as_one_connected_molecule() {
        let st = generate_structure("PC 16:0_18:1(9)").expect("should resolve");
        assert!(
            !st.regio_resolved,
            "sn assignment should be flagged as arbitrary"
        );
        assert!(
            !st.smiles.contains('.'),
            "should be one component: {}",
            st.smiles
        );
        assert!(
            !st.smiles.contains('*'),
            "should have no wildcards: {}",
            st.smiles
        );
        assert_chain_map_sane("PC 16:0_18:1(9)", &[(1, 16), (2, 18)]);

        // The stored form carries the same atoms plus the ambiguity.
        let plain = generate_smiles("PC 16:0_18:1(9)").expect("should resolve");
        assert!(
            plain.starts_with(&st.smiles) && plain.contains("swappable(sn1,sn2)"),
            "stored form should be the depiction plus the sn statement: {plain}"
        );
    }

    /// C1 is the carboxyl carbon, so the localized double bond declared at
    /// C9 must sit between the 9th and 10th mapped atoms.
    #[test]
    fn structure_map_carbon_numbering_lines_up_with_declared_db_position() {
        let st = generate_structure("FA 18:1(9)").expect("should resolve");
        let chain = &st.chains[0];
        let c9 = chain.carbons[8];
        let c10 = chain.carbons[9];
        let symbols = atom_symbols(&st.smiles);
        assert_eq!(symbols[c9], "C");
        assert_eq!(symbols[c10], "C");

        // C1 carries the carbonyl branch, so its token spans two atoms.
        assert_eq!(
            chain.carbons[1] - chain.carbons[0],
            2,
            "C1 should own its =O atom"
        );

        // Walk the string to the atom right after C9 and confirm the '='.
        let mut atom_idx = 0usize;
        let mut db_between_c9_c10 = false;
        let chars: Vec<char> = st.smiles.chars().collect();
        for (i, &c) in chars.iter().enumerate() {
            if matches!(c, 'C' | 'N' | 'O' | 'P' | 'S' | 'c' | 'n' | 'o' | '*') {
                if atom_idx == c9 && chars.get(i + 1) == Some(&'=') {
                    db_between_c9_c10 = true;
                }
                atom_idx += 1;
            }
        }
        assert!(db_between_c9_c10, "expected C9=C10 in {}", st.smiles);
    }

    #[test]
    fn tg_with_oxo_and_localized_db() {
        let s = generate_smiles("TG 18:0/18:0/18:1(9);5oxo").expect("should resolve");
        assert_balanced(&s);
        assert!(s.contains("C(=O)CCCC(=O)CCCC=CCCCCCCCC"));
        // Every position is known, so nothing is left to annotate.
        assert!(!s.contains(" |"), "{s}");
    }

    #[test]
    fn pc_with_geometry() {
        let s = generate_smiles("PC 16:0/18:1(9Z)").expect("should resolve");
        assert_balanced(&s);
        assert!(s.contains("[N+](C)(C)C"));
        assert!(s.contains("/C=C\\"));
    }

    #[test]
    fn pe_trans_geometry() {
        let s = generate_smiles("PE 16:0/18:1(9E)").expect("should resolve");
        assert_balanced(&s);
        assert!(s.contains("/C=C/"));
    }

    #[test]
    fn fa_unlocalized_db_gets_sg_blocks() {
        let s = generate_smiles("FA 20:4").expect("should resolve via Sg: blocks");
        assert_balanced(&s);
        assert!(
            s.contains(" |Sg:"),
            "unlocalized double bonds should be flagged with Sg:"
        );
        let sg_count = s.matches("Sg:n:").count();
        assert_eq!(
            sg_count, 5,
            "20:4 needs 4 unlocalized double bonds -> 5 Sg: markers (N+1)"
        );
        assert!(
            s.contains("a+b+c+d+e="),
            "should have a 5-variable size constraint"
        );

        let expanded = smiles_expand(&s);
        assert!(!expanded.contains('|'));
        assert_balanced(&expanded);
        assert_eq!(
            expanded.chars().filter(|&c| c == 'C').count(),
            20,
            "expanded chain should have all 20 carbons"
        );
        assert_eq!(
            expanded.matches("C=C").count(),
            4,
            "all 4 double bonds should survive expansion"
        );
    }

    #[test]
    fn fa_fully_saturated_resolves() {
        let s = generate_smiles("FA 16:0").expect("should resolve");
        assert_balanced(&s);
        assert_eq!(s, "OC(=O)CCCCCCCCCCCCCCC");
    }

    #[test]
    fn amp_fa_hete() {
        let s = generate_smiles("AMP-FA 20:4(5,8,11,14);15OH").expect("should resolve");
        assert_balanced(&s);
        assert!(s.contains("[n+]1ccccc1"));
        assert!(s.contains("C(O)")); // localized hydroxyl
                                     // All 4 double-bond positions are known and their cis/trans geometry
                                     // is not, which a plain `C=C` already says: an undecorated double
                                     // bond is unspecified geometry in SMILES, so there is nothing to add.
        assert_eq!(s.matches("C=C").count(), 4);
        assert!(
            !s.contains(" |"),
            "fully positioned chain needs no blocks: {s}"
        );
    }

    #[test]
    fn shorthand2020_functional_group_order_is_accepted() {
        let s = generate_smiles("FA 20:4(5,8,11,14);15OH")
            .expect("Shorthand2020 hydroxyl syntax should resolve");
        assert_balanced(&s);
        assert!(s.contains("C(O)"));
    }

    /// The bracketed tail must be stripped before parsing, so a name with
    /// one behaves exactly as the same name without one — including when
    /// that answer is `None`.
    #[test]
    fn confidence_display_tail_is_ignored_by_structure_parser() {
        for tail in ["", " [DB sn1: Δ5 100%, Δ8 100% | Δ14 50%]"] {
            let s = generate_smiles(&format!("FA 20:4(5,8,11,14);15OH{tail}"))
                .expect("display tail must not make the canonical structure unparsable");
            assert_balanced(&s);
        }
        // The space remaining after the hydroxyl cannot hold four
        // unlocalized double bonds without exceeding the named chain length.
        assert_eq!(generate_smiles("FA 20:4;11OH"), None);
        assert_eq!(
            generate_smiles("FA 20:4;11OH [DB sn1: Δ5 100%]"),
            None,
            "tail-stripping must not change the answer"
        );
    }

    #[test]
    fn amp_fa_unlocalized_oxygen_returns_none() {
        // generic ;O with no OH/oxo/COOH breakdown -> still ambiguous, no
        // placeholder convention was requested for oxygen sites.
        assert!(generate_smiles("AMP-FA 20:4(5,8,11,14);O").is_none());
    }

    #[test]
    fn lpc_single_chain_free_oh() {
        let s = generate_smiles("LPC 16:0").expect("should resolve");
        assert_balanced(&s);
        assert!(s.contains("(O)"));
        assert!(!s.contains('|'));
    }

    #[test]
    fn cer_d18_1_16_0() {
        let s = generate_smiles("Cer d18:1(4)/16:0").expect("should resolve");
        assert_balanced(&s);
        assert!(s.contains("NC(=O)CCCCCCCCCCCCCCC"));
    }

    #[test]
    fn sm_d18_1_16_0() {
        let s = generate_smiles("SM d18:1(4)/16:0").expect("should resolve");
        assert_balanced(&s);
        assert!(s.contains("[N+](C)(C)C"));
    }

    #[test]
    fn slot_fragment_unlocalized_db_expands_without_corrupting_carbonyl_branch() {
        // The glycerol-ester prefix must shift the chain's Sg atom indexes,
        // keeping expansion outside the C1 carbonyl branch.
        for name in ["LPC 18:2", "PC 18:2/14:1", "PE 18:2/14:1", "LPE 18:2"] {
            let s = generate_smiles(name).expect("should resolve");
            assert_balanced(&s);
            let expanded = smiles_expand(&s);
            assert!(!expanded.contains('|'));
            assert_balanced(&expanded);
            assert!(
                !expanded.contains("=OC"),
                "{name}: padding must never land inside a (=O) branch: {expanded}"
            );
        }

        // LPC 18:2: glycerophosphocholine backbone (11 headgroup atoms +
        // 3 glycerol carbons + 1 free OH) + an 18-carbon acyl chain.
        let s = generate_smiles("LPC 18:2").unwrap();
        let expanded = smiles_expand(&s);
        assert_eq!(
            expanded.chars().filter(|&c| c == 'C').count(),
            3 + 5 + 18,
            "3 glycerol C + 5 phosphocholine C + 18 chain C"
        );
        assert_eq!(expanded.matches("C=C").count(), 2);
    }

    #[test]
    fn ce_18_1() {
        let s = generate_smiles("CE 18:1(9Z)").expect("should resolve");
        assert_balanced(&s);
        assert!(s.contains("/C=C\\"));
    }

    #[test]
    fn car_ester_oxygen_bonds_to_the_carbonyl_carbon() {
        // Acylcarnitine: R-C(=O)-O-CH(CH2COO-)-CH2-N+(CH3)3. The ester O
        // bonds directly to the chain's carbonyl carbon (C1).
        let s = generate_smiles("CAR 18:0").expect("should resolve");
        assert_balanced(&s);
        assert_eq!(s, "O(C(=O)CCCCCCCCCCCCCCCCC)[C@H](CC(=O)[O-])C[N+](C)(C)C");
        assert!(
            s.starts_with("O(C(=O)"),
            "ester O must be directly bonded to the carbonyl carbon"
        );

        let s2 = generate_smiles("CAR 18:1").expect("should resolve");
        assert_balanced(&s2);
        assert!(s2.starts_with("O(C(=O)"));
        let expanded = smiles_expand(&s2);
        assert_balanced(&expanded);
        assert!(expanded.starts_with("O(C(=O)"));
        // 18 chain carbons + carnitine's own 7 ([C@H] + CH2-C(=O)O- branch (2) + CH2 + N(CH3)3 (3)).
        assert_eq!(expanded.chars().filter(|&c| c == 'C').count(), 18 + 7);
    }

    #[test]
    fn cl_unlocalized_db_gets_sg_regions() {
        // LipidOracle's EAD engines always join CL/DG/TG chains with "_",
        // never "/".
        let s = generate_smiles("CL 18:2_18:2_18:2_18:2").expect("should resolve");
        assert_balanced(&s);
        assert!(
            s.contains("Sg:n:"),
            "unlocalized double bonds should be flagged with Sg:"
        );
        assert!(s.contains("swappable(sn1,sn2,sn3,sn4)"), "{s}");
    }

    #[test]
    fn cl_four_chains_saturated_still_unresolved_regiochemistry() {
        let s = generate_smiles("CL 18:0_18:0_18:0_18:0").expect("should resolve");
        assert_balanced(&s);
        // CL is always "_"-joined by this codebase's own convention, so
        // regiochemistry is unresolved even though every chain here is
        // fully saturated.
        assert!(s.ends_with("swappable(sn1,sn2,sn3,sn4)"), "{s}");
        for sn in ["sn1", "sn2", "sn3", "sn4"] {
            assert_eq!(
                s.matches(sn).count(),
                2,
                "{sn} labels an atom and is named: {s}"
            );
        }
        assert!(!s.contains("Sg:"));
    }

    #[test]
    fn pc_slash_joined_stays_resolved_even_if_incomplete() {
        // "/" always means resolved, regardless of chain count.
        let s = generate_smiles("PC 16:0/18:1(9Z)").expect("should resolve");
        assert!(!s.contains("swappable("));
    }

    #[test]
    fn pc_underscore_joined_with_unlocalized_chains_keeps_all_ambiguity() {
        let s = generate_smiles("PC 16:1_18:1").expect("should resolve via Sg: blocks");
        assert_balanced(&s);
        assert!(
            s.contains("Sg:n:"),
            "unlocalized double bonds should be flagged with Sg:"
        );
        assert!(s.contains("swappable(sn1,sn2)"), "{s}");
        assert_eq!(
            s.matches('*').count(),
            0,
            "no R-group sites, so no wildcards: {s}"
        );
    }

    #[test]
    fn unknown_headgroup_returns_none() {
        assert!(generate_smiles("PIP2 16:0/18:1(9)").is_none());
    }

    #[test]
    fn placeholder_positions_avoid_existing_and_neighbors() {
        let extra = placeholder_positions(20, &[9], 1);
        // 9 and its neighbor 10 are blocked; 11 is the next free slot.
        assert_eq!(extra, vec![11]);
        assert!(!extra.contains(&9) && !extra.contains(&10));
        let extra2 = placeholder_positions(8, &[], 2);
        // short chain: must still avoid adjacency between the two picks
        assert_eq!(extra2.len(), 2);
        assert!((extra2[0] as i64 - extra2[1] as i64).abs() >= 2);
    }

    #[test]
    fn class_needs_multi_chain_matches_known_multi_and_single_chain_classes() {
        // LipidOracle's EAD engines rely on this exact list
        // to reject sum-composition parents (e.g. "PC O-32:1", no chain
        // split) before attempting to localize a double bond "within" the
        // sum as if it were one real chain -- keep this in sync with
        // generate_smiles's own shorthand rejection above.
        for class in [
            "DG", "TG", "CL", "PC", "PE", "PS", "PG", "PI", "PA", "Cer", "CerP", "SM", "HexCer",
            "IPC",
        ] {
            assert!(
                class_needs_multi_chain(class),
                "{class} should need multiple chains"
            );
        }
        for class in [
            "FA", "CE", "CAR", "NAE", "LPC", "LPE", "LPS", "LPG", "LPI", "LPA", "ST", "S1P", "Sph",
            "SB",
        ] {
            assert!(
                !class_needs_multi_chain(class),
                "{class} should not need multiple chains"
            );
        }
    }

    #[test]
    fn shorthand_pc_19_2_rejected() {
        // "PC 19:2" is ambiguous: total composition without explicit chains.
        // PC requires 2 chains but gets only 1 shorthand token "19:2",
        // so this should return None instead of creating a fake 19-carbon chain.
        assert!(generate_smiles("PC 19:2").is_none());
    }

    #[test]
    fn shorthand_tg_54_3_rejected() {
        // "TG 54:3" is ambiguous: total composition for 3 chains.
        // Should return None.
        assert!(generate_smiles("TG 54:3").is_none());
    }

    #[test]
    fn shorthand_dg_36_2_rejected() {
        // "DG 36:2" is ambiguous: total composition for 2 chains.
        // Should return None.
        assert!(generate_smiles("DG 36:2").is_none());
    }

    #[test]
    fn explicit_pc_16_0_slash_18_1_accepted() {
        // "PC 16:0/18:1" has explicit chains with "/" separator, so it's valid.
        let s = generate_smiles("PC 16:0/18:1").expect("explicit chains should work");
        assert_balanced(&s);
        assert!(!s.contains("swappable("));
    }

    #[test]
    fn explicit_pc_16_1_underscore_18_1_accepted() {
        // "PC 16:1_18:1" has explicit chains with "_" separator, so it's valid.
        let s = generate_smiles("PC 16:1_18:1").expect("explicit chains should work");
        assert_balanced(&s);
        // Both chains are unlocalized, so this takes the Sg:-keeping path.
        assert!(s.contains("Sg:n:"), "{s}");
    }

    #[test]
    fn single_chain_fa_18_1_accepted() {
        // FA and other single-chain lipids should work with shorthand like "FA 18:1".
        let s = generate_smiles("FA 18:1").expect("single-chain shorthand should work");
        assert_balanced(&s);
    }

    #[test]
    fn fa_with_hydroxyl_at_known_position() {
        // FA with hydroxyl at position 15
        let s = generate_smiles("FA 20:4(5,8,11,14);15OH").expect("FA with OH should work");
        assert_balanced(&s);
        assert!(s.contains("C(O)"), "hydroxyl group should be present");
        // Everything about this chain is positioned, so it comes out as
        // plain SMILES: the four double bonds' geometry is unspecified,
        // which a bare `C=C` already conveys.
        assert_eq!(s, "OC(=O)CCCC=CCC=CCC=CCC=C(O)CCCCC");
    }

    #[test]
    fn fa_with_hydroxyl_at_multiple_positions() {
        // FA with hydroxyls at positions 3 and 5
        let s = generate_smiles("FA 18:1(9);3OH,5OH").expect("FA with multiple OH should work");
        assert_balanced(&s);
        let oh_count = s.matches("C(O)").count();
        assert_eq!(oh_count, 2, "should have 2 hydroxyl groups");
    }

    #[test]
    fn fa_with_ketone_at_known_position() {
        // FA with ketone (oxo) at position 5
        let s = generate_smiles("FA 18:1(9);5oxo").expect("FA with oxo should work");
        assert_balanced(&s);
        assert!(s.contains("C(=O)"), "ketone group should be present");
    }

    #[test]
    fn fa_with_hydroxyls_and_ketones() {
        // FA with both hydroxyls and ketones
        let s = generate_smiles("FA 20:2(5,11);3OH;8oxo").expect("FA with OH and oxo should work");
        assert_balanced(&s);
        assert!(s.contains("C(O)"), "hydroxyl should be present");
        assert!(s.contains("C(=O)"), "ketone should be present");
    }

    #[test]
    fn tg_with_modifications_on_one_chain() {
        // TG with one chain having modifications
        let s = generate_smiles("TG 18:0/18:1(9);3OH/18:0")
            .expect("TG with chain modification should work");
        assert_balanced(&s);
        assert!(s.contains("C(O)"), "chain with hydroxyl should be present");
    }

    #[test]
    fn dg_with_hydroxyl_unresolved_regiochemistry() {
        // DG with hydroxyl on unresolved chains
        let s =
            generate_smiles("DG 18:1(9);3OH_18:0").expect("DG with OH and unresolved should work");
        assert_balanced(&s);
        assert!(
            s.contains("swappable(sn1,sn2)"),
            "unresolved regiochemistry needs a swappable token: {s}"
        );
        assert!(s.contains("C(O)"), "hydroxyl should be present");
    }

    #[test]
    fn generic_oxygen_without_position_returns_none() {
        // A generic oxygen count has no structural group or candidate site.
        assert!(
            generate_smiles("FA 18:1(9);O").is_none(),
            "generic O without position should be rejected"
        );
        assert!(
            generate_smiles("FA 18:1(9);O2").is_none(),
            "generic O2 without position should be rejected"
        );
    }

    #[test]
    fn hydroxyl_with_unlocalized_double_bond() {
        // A hydroxyl far enough along the chain to leave room for a double
        // bond on either side of it. The name does not say how the two bonds
        // distribute across the split and no CXSMILES this generator writes
        // can leave that open, so the only honest answer is to refuse.
        assert!(
            generate_smiles("FA 18:2;5OH").is_none(),
            "a localized group that splits the chain must be rejected"
        );

        // The same group on C3 pins nothing: there is no room for a double
        // bond before it, so the flexible run after it is the whole story.
        let s = generate_smiles("FA 18:2;3OH")
            .expect("a group with no room in front of it should still resolve");
        assert_balanced(&s);
        assert!(s.contains("C(O)"), "hydroxyl should be present");
        assert!(
            s.contains(" |Sg:"),
            "unlocalized double bonds should be flagged with Sg:"
        );

        let expanded = smiles_expand(&s);
        assert!(!expanded.contains('|'));
        assert_balanced(&expanded);
        assert!(
            expanded.contains("C(O)"),
            "hydroxyl should survive expansion"
        );
        assert_eq!(expanded.chars().filter(|&c| c == 'C').count(), 18);
        assert_eq!(expanded.matches("C=C").count(), 2);
    }

    /// Every Table 1A group this generator renders must actually reach the
    /// string. A group silently dropped still produces a plausible-looking
    /// SMILES — the whole chain, minus the thing the name was about — which
    /// is why this asserts the branch is present rather than just that the
    /// name resolved.
    #[test]
    fn every_table_1a_substituent_reaches_the_smiles() {
        for (abbr, branch) in SUBSTITUENTS {
            let name = format!("FA 18:0;5{abbr}");
            let s = generate_smiles(&name)
                .unwrap_or_else(|| panic!("{name} should resolve: Table 1A group {abbr}"));
            assert_balanced(&s);
            assert!(
                s.contains(&format!("C({branch})")),
                "{name}: {abbr} should appear as C({branch}) in {s}"
            );
            assert_chain_map_sane(&name, &[(1, 18)]);
        }
    }

    /// The Table 1A groups whose shorthand names only half a structure. Each
    /// has to be refused outright: emitting the chain without them would
    /// lose the modification silently, and inventing the other half is the
    /// fabrication this crate exists to avoid.
    #[test]
    fn underspecified_table_1a_groups_are_refused() {
        for abbr in UNRENDERABLE {
            let name = format!("FA 18:0;5{abbr}");
            assert_eq!(generate_smiles(&name), None, "{name} should be refused");
        }
        // C1 of an acyl chain is the ester carbonyl; `1OMe` is the methyl
        // ester, a different linkage rather than a substituted carbon.
        assert_eq!(generate_smiles("FA 18:0;1OMe"), None);
        // An unknown abbreviation is a typo, not a licence to ignore it.
        assert_eq!(generate_smiles("FA 18:0;5Zz"), None);
    }

    /// Table 1B rings, in the paper's `[start-end cyX:DBE(pos)]` spelling,
    /// with the three named acids it gives as examples.
    #[test]
    fn table_1b_rings_close_across_the_named_carbons() {
        // Lactobacillic acid: cyclopropane spanning C11-C13.
        let lacto = generate_smiles("FA 19:0;[11-13cy3:0]").expect("should resolve");
        assert_balanced(&lacto);
        assert_eq!(
            lacto.matches("%10").count(),
            2,
            "one open, one close: {lacto}"
        );
        assert!(lacto.contains("CC%10CC%10C"), "{lacto}");

        // Sterculic acid: the ring carries a double bond of its own, which
        // is additional to the chain's `19:0` count.
        let sterculic = generate_smiles("FA 19:0;[9-11cy3:1(9)]").expect("should resolve");
        assert!(sterculic.contains("C%10=CC%10"), "{sterculic}");

        // Gorlic acid: a cyclopentene ring plus a chain double bond.
        let gorlic = generate_smiles("FA 18:1(6Z);[14-18cy5:1(15)]").expect("should resolve");
        assert!(gorlic.contains("/C=C\\"), "chain DB kept: {gorlic}");
        assert!(gorlic.contains("C%10C=CCC%10"), "{gorlic}");

        for name in [
            "FA 19:0;[11-13cy3:0]",
            "FA 19:0;[9-11cy3:1(9)]",
            "FA 18:1(6Z);[14-18cy5:1(15)]",
        ] {
            assert_chain_map_sane(name, &[(1, if name.contains("18:1") { 18 } else { 19 })]);
        }

        // A span that disagrees with the ring size, or runs off the chain,
        // is malformed rather than something to round off.
        assert_eq!(generate_smiles("FA 19:0;[11-13cy5:0]"), None);
        assert_eq!(generate_smiles("FA 19:0;[18-20cy3:0]"), None);
    }

    /// An epoxide is a ring too: two adjacent carbons bridged by an oxygen,
    /// closed by the same `%nn` label. Both accepted epoxide spellings agree.
    #[test]
    fn epoxide_renders_as_a_three_membered_ring() {
        let s = generate_smiles("FA 18:0;5Ep").expect("should resolve");
        assert_balanced(&s);
        assert!(s.contains("C%10OC%10"), "{s}");
        assert_eq!(generate_smiles("FA 18:0;ep(5)").as_deref(), Some(&s[..]));
        assert_chain_map_sane("FA 18:0;5Ep", &[(1, 18)]);
    }

    /// Two ring-bearing chains in one molecule must not close each other's
    /// rings, and neither may collide with the fixed digits a headgroup
    /// template already uses (the sterol nucleus alone spends 1, 2 and 3).
    #[test]
    fn ring_labels_stay_distinct_across_chains_and_templates() {
        let tg = generate_smiles("TG 18:0;5Ep/18:0;9Ep/18:0").expect("should resolve");
        assert_balanced(&tg);
        assert_eq!(tg.matches("%10").count(), 2, "{tg}");
        assert_eq!(tg.matches("%11").count(), 2, "{tg}");
        assert_chain_map_sane("TG 18:0;5Ep/18:0;9Ep/18:0", &[(1, 18), (2, 18), (3, 18)]);

        // The sterol template's own 1/2/3 are untouched by the chain's %10.
        let ce = generate_smiles("CE 18:0;9Ep").expect("should resolve");
        assert_balanced(&ce);
        assert!(ce.contains("C12(") && ce.contains("C%10OC%10"), "{ce}");
    }

    /// The paper's combined example, which puts a ring, groups inside the
    /// ring's brackets, and a group outside it all on one chain.
    #[test]
    fn ring_brackets_may_carry_their_own_functional_groups() {
        let s = generate_smiles("FA 20:2(5,8);[11-15cy5;13OH];18OH").expect("should resolve");
        assert_balanced(&s);
        assert_eq!(s.matches("C(O)").count(), 2, "both hydroxyls: {s}");
        assert_eq!(s.matches("%10").count(), 2, "{s}");
        assert_chain_map_sane("FA 20:2(5,8);[11-15cy5;13OH];18OH", &[(1, 20)]);
    }

    /// `(OH)2` is the paper's form for "two of these, positions unknown".
    /// Each one needs its own `m:` block — collapsing them to one would
    /// under-report the group count.
    #[test]
    fn parenthesized_multiplicity_gives_one_m_block_per_group() {
        let s = generate_smiles("FA 20:3(5,8,11);(OH)2").expect("should resolve");
        assert_eq!(s.matches("m:").count(), 2, "{s}");
        assert_eq!(
            s.split(" |").next().unwrap().matches(".*O").count(),
            2,
            "{s}"
        );
    }

    /// Positions in front, comma-separated, several groups per name — the
    /// spelling the paper actually specifies.
    #[test]
    fn positions_precede_abbreviations_and_stack() {
        // Phytanic acid: four methyl branches on a 16-carbon chain.
        let phytanic = generate_smiles("FA 16:0;3Me,7Me,11Me,15Me").expect("should resolve");
        assert_eq!(phytanic.matches("C(C)").count(), 4, "{phytanic}");
        assert_chain_map_sane("FA 16:0;3Me,7Me,11Me,15Me", &[(1, 16)]);

        let mixed = generate_smiles("FA 20:3(5Z,13E,17E);11OH,15OH;9oxo").expect("should resolve");
        assert_eq!(mixed.matches("C(O)").count(), 2, "{mixed}");
        assert_eq!(
            mixed.matches("C(=O)").count(),
            2,
            "C1 plus the 9-oxo: {mixed}"
        );
    }

    #[test]
    fn simple_fa() {
        let s = generate_smiles("FA 18:1(9)").expect("should resolve");
        assert_balanced(&s);
        assert!(s.contains("C=C"), "double bond should be present");
    }

    #[test]
    fn expand_cxsmiles_resolves_to_full_chain_length() {
        // FA 18:2: 18 total carbons, 2 unlocalized double bonds.
        let cxsmiles = generate_smiles("FA 18:2").expect("should resolve");
        let expanded = smiles_expand(&cxsmiles);

        assert!(
            !expanded.contains('|'),
            "expanded SMILES should drop the CXSMILES suffix"
        );
        assert_balanced(&expanded);

        let carbon_count = expanded.chars().filter(|&c| c == 'C').count();
        assert_eq!(
            carbon_count, 18,
            "expanded chain should have all 18 carbons: {}",
            expanded
        );
        assert_eq!(
            expanded.matches("C=C").count(),
            2,
            "both double bonds should survive expansion: {}",
            expanded
        );
    }

    #[test]
    fn expand_cxsmiles_handles_multiple_unlocalized_dbs() {
        let cxsmiles = "CC=CC=CC(=O)O |Sg:n:0:a:ht,Sg:n:2:b:ht,Sg:n:4:c:ht| a+b+c=15";
        let expanded = smiles_expand(cxsmiles);
        assert!(!expanded.contains('|'));
        assert_balanced(&expanded);
        assert_eq!(expanded.chars().filter(|&c| c == 'C').count(), 18);
    }

    #[test]
    fn expand_cxsmiles_passthrough_without_sg_blocks() {
        // Plain SMILES with no CXSMILES suffix is returned unchanged.
        let s = generate_smiles("FA 16:0").expect("should resolve");
        assert_eq!(smiles_expand(&s), s);
    }

    #[test]
    fn fa_unlocalized_db_has_cx_fields() {
        // FA with unlocalized double bonds
        let s = generate_smiles("FA 18:2").expect("should resolve");
        let smiles_part = s.split(' ').next().unwrap_or(&s);
        assert_balanced(smiles_part);

        // Should use Sg notation format: OC(=O)CC=CC=CC |Sg:n:3:a:ht,Sg:n:5:b:ht,Sg:n:8:c:ht| a+b+c=N
        // (positions 3/5/8, not 0/2/4: the free acid's O and the carbonyl
        // C(=O) both precede the Sg-marked tail and count as atoms too; the
        // final marker is on the terminal carbon after the last C=C.)
        assert!(s.contains("C=C"), "should have explicit double bonds");
        assert!(
            smiles_part.ends_with("C=CC"),
            "uncertain unsaturated chains must end in C=CC"
        );
        assert!(
            s.contains("Sg:n:3:a:ht"),
            "should have first Sg after the O/C(=O) prefix"
        );
        assert!(
            s.contains("Sg:n:5:b:ht"),
            "should have second Sg two atoms later"
        );
        assert!(
            s.contains("Sg:n:8:c:ht"),
            "should have final Sg on the terminal carbon"
        );
        assert!(
            s.contains("|Sg:"),
            "should have CXSMILES suffix with Sg blocks"
        );
        assert!(
            s.starts_with("OC(=O)"),
            "should have the free acid's hydroxyl and carboxyl group"
        );

        assert!(!s.contains("[H]"), "should NOT have [H] prefix");

        // Verify constraint format: a+b+c+... = totalC - #DBs - 2
        assert!(s.contains("a+b+c="), "should have constraint a+b+c=");
        assert!(
            s.contains("a+b+c=14"),
            "FA 18:2 should have constraint a+b+c=14 (18-2-2)"
        );
    }

    #[test]
    fn single_unlocalized_db_has_terminal_sg_marker() {
        let s = generate_smiles("FA 18:1").expect("should resolve");
        let smiles_part = s.split(' ').next().unwrap_or(&s);

        assert!(smiles_part.ends_with("C=CC"));
        assert!(s.contains("Sg:n:3:a:ht,Sg:n:6:b:ht"));
        assert!(
            s.contains("a+b=15"),
            "the terminal carbon reduces the variable count by one"
        );

        let expanded = smiles_expand(&s);
        assert_eq!(expanded.chars().filter(|&c| c == 'C').count(), 18);
        assert_eq!(expanded.matches("C=C").count(), 1);
    }

    #[test]
    fn unknown_geometry_is_just_a_plain_double_bond() {
        let unspecified = generate_smiles("FA 18:1(9);OH").expect("should resolve");
        assert!(unspecified.contains("|m:"));
        assert!(unspecified.contains("C=C"), "{unspecified}");
        assert!(!unspecified.contains('/') && !unspecified.contains('\\'));

        let explicit = generate_smiles("FA 18:1(9Z);OH").expect("should resolve");
        assert!(explicit.contains("/C=C\\"));
    }

    #[test]
    fn pc_complex_ambiguities() {
        let test_name = "PC 16:2_18:1(7);5oxo";
        let s = generate_smiles(test_name).expect("should parse");

        assert_balanced(&s);
        assert!(s.contains("C(=O)C") || s.contains("(=O)C"), "{s}");
        assert!(
            s.contains("Sg:"),
            "should have Sg: for chain 1's unlocalized double bond"
        );
        assert!(s.contains("swappable(sn1,sn2)"), "{s}");
    }

    /// An unlocalized modification is a CXSMILES position-variation bond:
    /// `m:<dummy atom>:<candidate anchors>`. The anchors must be this
    /// chain's own carbons (minus the acyl C1), the variable end must be a
    /// `*` dummy carrying exactly one bond, and the component must not
    /// smuggle in a carbon the chain doesn't have.
    #[test]
    fn unlocalized_modification_lists_every_chain_carbon() {
        let s = generate_smiles("FA 18:1;OH").expect("should resolve");
        // OC(=O)CC=CC.*O -> chain carbons are atoms 1..=6, C1 (atom 1) is
        // the carboxyl, and atom 7 is the `*` the hydroxyl hangs off.
        assert!(s.contains("m:7:3.4.5.6"), "{s}");
        assert_eq!(s.split(" |").next().unwrap(), "OC(=O)CC=CC.*O", "{s}");

        // A ketone converts a carbon rather than adding one; an extra
        // carboxyl brings its own.
        let oxo = generate_smiles("FA 18:0;oxo").expect("should resolve");
        assert!(oxo.split(" |").next().unwrap().ends_with(".*=O"), "{oxo}");
        let cooh = generate_smiles("FA 18:0;COOH").expect("should resolve");
        assert!(
            cooh.split(" |").next().unwrap().ends_with(".*C(=O)O"),
            "{cooh}"
        );
    }

    /// The position-variation stub survives into the depiction form
    /// unchanged and adds no carbon the molecule does not have.
    #[test]
    fn depiction_keeps_the_position_variation_stub() {
        for (name, expected) in [
            ("FA 18:1;OH", "OC(=O)CCCCCCCC=CCCCCCCCC.*O"),
            ("FA 18:0;oxo", "OC(=O)CCCCCCCCCCCCCCCCC.*=O"),
            ("FA 18:0;COOH", "OC(=O)CCCCCCCCCCCCCCCCC.*C(=O)O"),
        ] {
            let depicted = smiles_expand(&generate_smiles(name).unwrap());
            assert_eq!(depicted, expected, "{name}");
        }
        // A charged bracket atom inside a real structure is not a
        // position-variation placeholder and must survive untouched.
        let pc = smiles_expand(&generate_smiles("PC 16:0/18:1(9)").unwrap());
        assert!(pc.contains("[O-]") && pc.contains("[N+]"), "{pc}");

        // The stub sits between the chain and whatever follows, so every
        // chain's carbon indices still have to line up around it.
        assert_chain_map_sane("FA 18:1;OH", &[(1, 18)]);
        assert_chain_map_sane("Cer d18:1(4)/16:0;OH", &[(1, 18), (2, 16)]);
        assert_chain_map_sane("DG 18:1(9);5OH_18:1;OH", &[(1, 18), (2, 18)]);
    }

    /// An sn-unresolved name must lose nothing to the encoding: both chains
    /// keep every carbon, and the ambiguity rides alongside rather than
    /// displacing them.
    #[test]
    fn swappable_names_keep_every_carbon() {
        let s = generate_smiles("DG 18:1(9);5OH_18:0").expect("should resolve");
        assert_balanced(&s);
        assert!(s.ends_with("swappable(sn1,sn2)"), "{s}");
        assert_chain_map_sane("DG 18:1(9);5OH_18:0", &[(1, 18), (2, 18)]);
    }

    /// A chain can be both unlocalized and swappable.
    #[test]
    fn sg_and_swappable_coexist() {
        let s = generate_smiles("PC 16:0_18:1").expect("should resolve");
        assert_balanced(&s);
        assert!(s.contains("Sg:n:"), "unlocalized double bond kept: {s}");
        assert!(
            s.contains("swappable(sn1,sn2)"),
            "sn ambiguity kept too: {s}"
        );

        // And the same for an unlocalized modification.
        let m = generate_smiles("PC 16:0_18:1(9);OH").expect("should resolve");
        assert!(m.contains("m:"), "{m}");
        assert!(m.contains("swappable(sn1,sn2)"), "{m}");
    }

    /// Every Sg-bearing chain expands to its declared carbon count.
    #[test]
    fn sg_chains_keep_every_carbon() {
        assert_chain_map_sane("DG 18:1(9);5OH_18:1;OH", &[(1, 18), (2, 18)]);
    }

    #[test]
    fn expand_cxsmiles_handles_multi_chain_offset_and_renaming() {
        // Two chains each with their own unlocalized double bonds, joined
        // with unresolved regiochemistry: exercises CxBuilder's atom
        // offsetting and variable renaming across multiple fragments, and
        // smiles_expand's positional (not name-based)
        // matching of Sg positions back to their constraint equation.
        let s = generate_smiles("PC 16:2_18:1").expect("should resolve");
        assert_balanced(&s);
        assert!(s.contains("swappable(sn1,sn2)"), "{s}");
        assert_eq!(
            s.matches("Sg:n:").count(),
            5,
            "chain1 (16:2, 2 unlocalized DBs -> 3 Sg) + chain2 (18:1, 1 unlocalized DB -> 2 Sg)"
        );

        let expanded = smiles_expand(&s);
        assert!(!expanded.contains('|'));
        assert_balanced(&expanded);
        // Backbone (glycerol C1-C3 + phosphocholine's 5 C's) plus every
        // carbon of both chains.
        assert_eq!(
            expanded.chars().filter(|&c| c == 'C').count(),
            3 + 5 + 16 + 18
        );
        assert_eq!(
            expanded.matches('*').count(),
            0,
            "chains are esterified in place, so no attachment stubs remain"
        );
        assert_eq!(
            expanded.matches("C=C").count(),
            3,
            "2 DBs on chain1 + 1 DB on chain2"
        );
    }

    #[test]
    fn pc_with_variable_db_multichain() {
        // PC with one chain having variable DBs and one known chain
        // PC 18:4(5,8,11,14);3oxo;16OH_18:2;oxo;OH
        // Chain 1: fully known (variable DBs + known modifications)
        // Chain 2: variable DBs + variable modifications
        let s = generate_smiles("PC 18:4(5,8,11,14);3oxo;16OH_18:2;oxo;OH")
            .expect("should resolve PC with variable content");
        let blocks = s.split(" |").nth(1).unwrap_or("");
        assert!(blocks.contains("Sg:n:"), "chain 2 has unlocalized DBs: {s}");
        assert!(
            blocks.contains("m:"),
            "chain 2 has unlocalized oxygens: {s}"
        );
        assert!(s.contains("swappable(sn1,sn2)"), "{s}");
        // The `*` present is the position-variation stub, not an attachment
        // point for a chain.
        assert!(s.contains(".*O"), "{s}");
    }

    #[test]
    fn extended_classes_support_variable_content() {
        for name in [
            "PE 18:2_18:1(9)",
            "PS 18:2_18:1(9)",
            "PG 18:2_18:1(9)",
            "PA 18:2_18:1(9)",
            "LPE 18:2",
            "LPS 18:2",
        ] {
            let cxsmiles = generate_smiles(name).unwrap_or_else(|| panic!("{name} should resolve"));
            assert!(cxsmiles.contains("Sg:"), "{name}: {cxsmiles}");
        }
    }
}
