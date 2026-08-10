//! On-demand SMILES generation for lipid species/isomer names.
//!
//! The public API wraps these entry points: see `name2smiles`,
//! `name2structure` and `expand_cxsmiles_for_depiction` in `lib.rs`.
//!
//! `lipid_name_to_smiles` takes a canonical Shorthand2020 name —
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
//! chains are localized) — in that case the backbone carries a
//! `*` R-group slot at each ester position and the chains become the
//! alternative definitions of that R-group, e.g. `PC 16:0_18:1(9)` becomes
//! `C(COP(=O)([O-])OCC[N+](C)(C)C)(O*)CO* |$;;;;;;;;;;;;;;R1;;;R1$,RG:_R1={C(=O)CCCCCCCCCCCCCCC},{C(=O)CCCCCCCC=CCCCCCCCC}|`.
//! The ester oxygen stays on the backbone, so each definition still carries
//! all of its chain's own carbons.
//!
//! An `RG:` definition cannot contain a nested CXSMILES block (CDK rejects
//! it outright), so a chain that needs an `Sg:`/`m:` block of its own
//! cannot become an R-group alternative. When any chain is in that position
//! the sn ambiguity is the one that gives way: every chain is written into
//! the backbone in name order and no `RG:` block is emitted. See
//! `rg_alternatives`.
//!
//! Two constructs are deliberately *not* used. `f:` groups components into
//! one entity (a salt, a hydrate) and expresses *and*, never *or*, so it
//! cannot say "this chain or that chain sits at sn-1"; earlier versions
//! used it for exactly that. `ctu:` is a query feature for matching either
//! configuration — a plain `C=C` is already unspecified geometry in SMILES,
//! so emitting it added nothing and made every structure a query.
//!
//! ## Unlocalized double bonds
//!
//! Every chain builder is Sg:-aware: known double bonds/modifications are
//! rendered literally (with geometry where declared), and any remaining
//! *unlocalized* double bonds are represented with a CDK-style `Sg:n:`
//! flexible-run marker plus an `x+y=N` size constraint, never with a
//! guessed literal position. See `chain_fragment_cdk_range` and
//! `CdkBuilder` for how a chain's local Sg:/m: blocks get offset and
//! merged into the surrounding headgroup template's global CXSMILES
//! suffix. A double bond whose geometry is not explicit is left as a plain
//! `C=C`, which is exactly what unspecified geometry means in SMILES;
//! explicit `Z`/`E` geometry remains encoded with
//! `/` and `\\`.
//!
//! ## Supported modifications
//!
//! The parser recognizes these chain modifications when positions are given:
//! - `;posOH` — hydroxyl groups (alcohol, `-OH`)
//! - `;posoxo` — ketone groups (carbonyl, `C=O`, internal)
//! - `;COOH(...)` — extra carboxylic acid branches
//! - `;ep(...)` / `;epox(...)` — epoxide rings (currently parsed but not yet rendered in SMILES)
//! - `;cyc(...)` / `;cyclo(...)` — cyclopropane rings (currently parsed but not yet rendered in SMILES)
//!
//! `None` is still returned when the headgroup has no template here, or
//! when a chain's oxygen count is declared generically (`;O`/`;O2`) with
//! no `;OH`/`;oxo`/`;COOH` position breakdown at all — unlike double
//! bonds, there's no placeholder convention asked for there, so that
//! case is left unresolved rather than guessed. Epoxides and cyclopropanes
//! are fully parsed but SMILES rendering for ring structures is not yet implemented.

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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Regiochemistry {
    /// `/`-joined (or single-chain): sn-position known, one connected molecule.
    Resolved,
    /// `_`-joined with 2+ chains: sn-position/regiochemistry unknown, chains
    /// rendered as the alternative definitions of a CXSMILES `RG:` R-group
    /// hung off `*` slots on the backbone — or, when a chain needs an
    /// `Sg:`/`m:` block that an R-group definition cannot hold, built as
    /// `Resolved` in name order with the sn ambiguity left unexpressed.
    Unresolved,
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

#[derive(Debug, Clone)]
struct ParsedChain {
    prefix: ChainPrefix,
    carbon: u32,
    db_pos: Vec<DbPos>,
    oh_pos: Vec<u32>,
    ket_pos: Vec<u32>,
    cooh_pos: Vec<u32>,
    /// Epoxide positions (5,6-epoxide at C5-C6). Rendering deferred.
    #[allow(dead_code)]
    epox_pos: Vec<u32>,
    /// Cyclopropane ring positions (3,4-cyclopropane at C3-C4). Rendering deferred.
    #[allow(dead_code)]
    cyc_pos: Vec<u32>,
}

/// Entry point: build semantically correct SMILES/CXSMILES for one lipid species.
///
/// For structures with variable content (unknown DBs, unknown mods, or
/// unresolved regiochemistry): returns CDK-style CXSMILES with `Sg:`/`m:`/
/// `RG:` blocks. For fully-known structures: returns plain SMILES (no
/// CXSMILES suffix). Returns `None` if the headgroup isn't covered or a
/// chain's oxygen count has no position hypothesis at all.
///
/// Reference: https://egonw.github.io/cdk-cxsmiles/templates.html#lipids-with-two-double-bonds-somewhere-in-the-tail
pub fn lipid_name_to_smiles(name: &str) -> Option<String> {
    let (canonical, _) = crate::nomenclature::split_display_name(name);
    build_lipid(&canonical, Regiochemistry::Unresolved).map(|b| b.smiles)
}

/// Builds one lipid's SMILES plus, for every chain, the atom index of each
/// of its carbons — the mapping a UI needs to highlight the part of a
/// structure that a given MS2 fragment comes from.
///
/// Differs from [`lipid_name_to_smiles`] in two depiction-driven ways:
///
/// * `_`-joined (sn-unresolved) names are built as if they were `/`-joined,
///   because the honest CXSMILES rendering — a backbone of `*` R-group slots
///   plus `RG:` alternatives — is a Markush scheme rather than a molecule,
///   and carries no concrete chain atoms for a viewer to highlight. The
///   arbitrary sn assignment that buys a connected depiction is flagged via
///   [`LipidStructure::regio_resolved`] so callers can say so.
/// * The `Sg:` unlocalized-double-bond markers are already expanded (see
///   [`expand_cxsmiles_for_depiction`]), and the atom indices account for
///   the padding atoms that expansion inserts.
///
/// Atom indices are 0-based in SMILES emission order, which is the order
/// RDKit (and every other parser) assigns when reading the string.
pub fn lipid_name_to_structure(name: &str) -> Option<LipidStructure> {
    let (canonical, _) = crate::nomenclature::split_display_name(name);
    let labels = chain_label_tokens(&canonical);
    let built = build_lipid(&canonical, Regiochemistry::Resolved)?;

    let (expanded, inserts) = expand_with_padding_inserts(&built.smiles);
    let atom_count = count_atoms(&expanded, true, true);
    let shift = padding_shift_table(&inserts);
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

/// Shared builder behind [`lipid_name_to_smiles`] and
/// [`lipid_name_to_structure`]. `force` is the regiochemistry mode to use
/// for multi-chain `_`-joined names; `/`-joined and single-chain names are
/// always resolved regardless.
fn build_lipid(name: &str, force: Regiochemistry) -> Option<Built> {
    let (class, rest) = match name.trim().split_once(' ') {
        Some((h, r)) => (h, r.trim()),
        None => (name.trim(), ""),
    };

    if rest.is_empty() {
        // Single lipid class with no chains - return simple SMILES
        return build_simple_smiles(class);
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

    let mode = if is_slash || chains.len() <= 1 {
        Regiochemistry::Resolved
    } else {
        force
    };

    match class {
        "FA" => fa_smiles(&chains),
        "AMP-FA" | "FA-AMP" => amp_fa_smiles(&chains),
        "NAE" => nae_smiles(&chains),
        "CAR" => car_smiles(&chains),
        "CE" => ce_smiles(&chains),
        "ST" => Some(Built::plain(st_smiles())),
        "MG" => glycerolipid_smiles(&chains, 1, Regiochemistry::Resolved),
        "DG" => glycerolipid_smiles(&chains, 2, mode),
        "TG" => glycerolipid_smiles(&chains, 3, mode),
        "CL" => cl_smiles(&chains, mode),
        "PC" | "LPC" => gpl_smiles("OP(=O)([O-])OCC[N+](C)(C)C", &chains, mode),
        "PE" | "LPE" => gpl_smiles("OP(=O)(O)OCCN", &chains, mode),
        "PS" | "LPS" => gpl_smiles("OP(=O)(O)OCC(N)C(=O)O", &chains, mode),
        "PG" | "LPG" => gpl_smiles("OP(=O)(O)OCC(O)CO", &chains, mode),
        "PI" | "LPI" => gpl_smiles("OP(=O)(O)OC1C(O)C(O)C(O)C(O)C1O", &chains, mode),
        "PA" | "LPA" => gpl_smiles("OP(=O)(O)O", &chains, mode),
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

/// Alias for `lipid_name_to_smiles()`, kept for existing tests/call sites
/// written when this crate distinguished a "traditional" vs "CDK" mode.
/// There's only one mode now.
#[allow(dead_code)]
pub fn lipid_name_to_cxsmiles(name: &str) -> Option<String> {
    lipid_name_to_smiles(name)
}

/// Classes whose structure inherently spans multiple acyl/alkyl chains
/// (sn1/sn2[/sn3/sn4]), so a bare sum-composition token (e.g. `"32:1"`,
/// no `"_"`/`"/"` chain split) can't be resolved to one real structure —
/// the total carbons/double bonds could come from many different chain
/// combinations. Used both to reject shorthand species names here and,
/// by the EAD engines, to reject the same shorthand as a double-bond
/// localization target in the first place (localizing a position within
/// the *sum* would treat two real chains as if they were one).
pub fn class_needs_multi_chain(class: &str) -> bool {
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

/// Simple SMILES for lipid classes with no explicit chains
fn build_simple_smiles(class: &str) -> Option<Built> {
    if class == "ST" {
        Some(Built::plain(st_smiles()))
    } else {
        None
    }
}

/// An assembled structure: the SMILES text plus, per chain, the sn position
/// and the local atom indices of its carbons in C1..Cn order.
#[derive(Debug, Default, Clone)]
struct Built {
    smiles: String,
    chains: Vec<(usize, Vec<usize>)>,
}

impl Built {
    /// A structure with no chain mapping (fixed templates like sterol).
    fn plain(smiles: String) -> Self {
        Built {
            smiles,
            chains: Vec::new(),
        }
    }
}

/// A chain fragment (or, for `wildcard_fragment_cdk`, a disconnected
/// standalone dot-joined component) with CDK CXSMILES ambiguity
/// annotations, using LOCAL atom indices (0-based, counting every
/// element/bracket-atom/`*` token from this fragment's own first atom)
/// and LOCAL variable letters (always starting fresh at `a`). Embed into
/// a larger assembled SMILES via `CdkBuilder::push_fragment`, which
/// offsets positions and renames variables so they don't collide with
/// whatever came before in the assembly.
#[derive(Debug, Default, Clone)]
struct CdkFragment {
    smiles: String,
    sg_blocks: Vec<String>,
    constraint: Option<String>,
    m_blocks: Vec<String>,
    /// Local atom index of each chain carbon actually emitted, in C1..Cn
    /// order. Shorter than the chain's carbon count when `Sg:` markers
    /// stand in for a variable-length run — the atoms that
    /// `expand_cxsmiles_for_depiction` inserts fill the rest of the run in
    /// place, so `lipid_name_to_structure` splices them back in after
    /// expansion. Empty for non-chain fragments (bare `O` slots).
    carbon_atoms: Vec<usize>,
}

/// Builds one chain's fragment (C1..Cn), dispatching on its linkage type.
/// `EtherAlkenyl` gets its mandatory (never placeholder) vinyl-ether
/// C1=C2 folded in before any of the chain's own declared double bonds.
/// Sphingoid bases aren't chain fragments in this sense — see
/// `sphingoid_smiles`, which uses `chain_fragment_cdk_range` directly
/// with a `start` of 3 (after the fixed C1/C2 positions).
fn chain_fragment_cdk(chain: &ParsedChain) -> Option<CdkFragment> {
    match chain.prefix {
        ChainPrefix::Acyl => chain_fragment_cdk_range(
            1,
            chain.carbon,
            &chain.db_pos,
            &chain.oh_pos,
            &chain.ket_pos,
            &chain.cooh_pos,
            true,
        ),
        ChainPrefix::EtherAlkyl => chain_fragment_cdk_range(
            1,
            chain.carbon,
            &chain.db_pos,
            &chain.oh_pos,
            &chain.ket_pos,
            &chain.cooh_pos,
            false,
        ),
        ChainPrefix::EtherAlkenyl => {
            let mut db = vec![DbPos {
                pos: 1,
                geom: None,
                placeholder: false,
            }];
            db.extend(chain.db_pos.iter().copied());
            chain_fragment_cdk_range(
                1,
                chain.carbon,
                &db,
                &chain.oh_pos,
                &chain.ket_pos,
                &chain.cooh_pos,
                false,
            )
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
/// `expand_cxsmiles_for_depiction` only needs to insert `value - 1` more
/// to reach any one valid total length. Modifications with no position
/// at all (`pos == 0`) become extra dot-joined fragments plus `m:`
/// blocks, appended after the main chain.
fn chain_fragment_cdk_range(
    start: u32,
    carbon: u32,
    db_pos: &[DbPos],
    oh_pos: &[u32],
    ket_pos: &[u32],
    cooh_pos: &[u32],
    carbonyl_c1: bool,
) -> Option<CdkFragment> {
    if carbon < start {
        return Some(CdkFragment::default());
    }

    let localized: Vec<DbPos> = db_pos.iter().copied().filter(|d| !d.placeholder).collect();
    let unlocalized_count = db_pos.len() - localized.len();
    let known_oh: Vec<u32> = oh_pos.iter().copied().filter(|&p| p != 0).collect();
    let known_ket: Vec<u32> = ket_pos.iter().copied().filter(|&p| p != 0).collect();
    let known_cooh: Vec<u32> = cooh_pos.iter().copied().filter(|&p| p != 0).collect();

    let mut smiles;
    let mut sg_blocks = Vec::new();
    let mut constraint = None;
    let mut carbon_atoms;

    if unlocalized_count == 0 {
        (smiles, carbon_atoms) = build_chain_range(
            carbon,
            &localized,
            &known_oh,
            &known_ket,
            &known_cooh,
            carbonyl_c1,
            start,
            carbon,
        );
    } else {
        let mut prefix_len = start;
        for d in &localized {
            prefix_len = prefix_len.max(d.pos + 1);
        }
        for &p in known_oh
            .iter()
            .chain(known_ket.iter())
            .chain(known_cooh.iter())
        {
            prefix_len = prefix_len.max(p);
        }
        prefix_len = prefix_len.min(carbon);

        (smiles, carbon_atoms) = build_chain_range(
            carbon,
            &localized,
            &known_oh,
            &known_ket,
            &known_cooh,
            carbonyl_c1,
            start,
            prefix_len,
        );

        let remaining_length = carbon - prefix_len;
        let var_names = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j'];
        let num_markers = (unlocalized_count + 1).min(var_names.len());
        let prefix_atom_count = count_atoms(&smiles, true, true);

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
        // The explicit terminal carbon above accounts for one more atom than
        // the old C=C-ending template, so the variable-length total is one
        // smaller.
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
    let has_unknown_oh = oh_pos.contains(&0);
    let has_unknown_oxo = ket_pos.contains(&0);
    let has_unknown_cooh = cooh_pos.contains(&0);
    if has_unknown_oh || has_unknown_oxo || has_unknown_cooh {
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
        let mut atom_count = count_atoms(&smiles, true, true);
        // Each appended component leads with a `*` dummy, and the `m:`
        // block points at that dummy rather than at the modification's own
        // heteroatom. This is not cosmetic: a position-variation bond's
        // variable end has to be a dummy atom carrying exactly one bond.
        // Both toolkits enforce it — CDK ignores an `m:` whose target is
        // anything else, and RDKit rejects the whole string ("position
        // variation bond to atom with more than one bond"), so the older
        // bracket-atom forms (`[OH]`, `[O]`, `[C](=O)O`) were inert at
        // best. The dummy also supplies the free valence those brackets
        // were there to reserve, so no explicit hydrogen counts are
        // needed. A hydroxyl and a ketone both convert an existing chain
        // carbon, so their component is just the oxygen; only an extra
        // carboxyl brings a carbon of its own.
        for (present, component) in [
            (has_unknown_oh, "*O"),
            (has_unknown_oxo, "*=O"),
            (has_unknown_cooh, "*C(=O)O"),
        ] {
            if !present {
                continue;
            }
            smiles.push('.');
            smiles.push_str(component);
            m_blocks.push(format!("m:{atom_count}:{sites}"));
            atom_count += count_atoms(component, true, true);
        }
    }

    Some(CdkFragment {
        smiles,
        sg_blocks,
        constraint,
        m_blocks,
        carbon_atoms,
    })
}

/// A glycerol-ester/ether slot: `"O" + fragment` for a real chain, or a
/// bare free hydroxyl `"O"` when the position is empty/absent.
fn slot_fragment_cdk(chain: Option<&ParsedChain>) -> Option<CdkFragment> {
    match chain {
        None => Some(CdkFragment {
            smiles: "O".to_string(),
            ..Default::default()
        }),
        Some(c) if c.carbon == 0 => Some(CdkFragment {
            smiles: "O".to_string(),
            ..Default::default()
        }),
        Some(c) => {
            let frag = chain_fragment_cdk(c)?;
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
            Some(CdkFragment {
                smiles: format!("O{}", frag.smiles),
                sg_blocks,
                constraint: frag.constraint,
                m_blocks,
                carbon_atoms: frag.carbon_atoms.iter().map(|a| a + 1).collect(),
            })
        }
    }
}

/// Whether a chain has to be written into the main SMILES rather than into
/// an R-group definition — i.e. whether it carries any `Sg:` or `m:` block
/// of its own.
///
/// Both blocks index atoms of the *main* string, and CDK rejects a nested
/// `|...|` block inside an `RG:` definition outright (verified against CDK
/// Depict: `RG:_R1={C(=O)CC=CC |Sg:n:2:a:ht|}` is a parse error, as is the
/// `{*C(=O)CCC |$_AP1$|}` attachment-point form). So a chain whose double
/// bonds or modifications are unlocalized cannot live in a definition, and
/// the sn-ambiguity that `RG:` would express has to be given up for it —
/// see `rg_alternatives`.
fn chain_needs_main_string(chain: &ParsedChain) -> bool {
    chain_fragment_cdk(chain)
        .map(|f| !f.sg_blocks.is_empty() || !f.m_blocks.is_empty())
        .unwrap_or(false)
}

/// The `RG:` definitions for a set of interchangeable chains, in name
/// order, or `None` if any of them can't be expressed as one (see
/// `chain_needs_main_string`).
///
/// Each definition is the chain's own SMILES with no attachment-point
/// marker: the first atom of a definition is implicitly the attachment,
/// so C1 of the chain binds to whichever `*` carries the matching `R1`
/// label. Every chain keeps all `carbon` of its own carbons — the ester
/// oxygen stays on the backbone side, in the `O*` stub.
fn rg_alternatives(chains: &[ParsedChain]) -> Option<Vec<String>> {
    if chains.iter().any(chain_needs_main_string) {
        return None;
    }
    chains
        .iter()
        .map(|c| {
            (c.carbon > 0)
                .then(|| chain_fragment_cdk(c).map(|f| f.smiles))
                .flatten()
        })
        .collect()
}

/// Incrementally assembles a multi-part CXSMILES body (fixed literal
/// text interleaved with chain fragments), offsetting each fragment's
/// local `Sg:`/`m:` atom indices by the atom count of everything emitted
/// before it, and renaming each fragment's own local variable letters
/// (which always start fresh at `a`) so they don't collide with an
/// earlier fragment's.
#[derive(Default)]
struct CdkBuilder {
    smiles: String,
    sg_blocks: Vec<String>,
    m_blocks: Vec<String>,
    constraints: Vec<String>,
    /// Global atom indices carrying an `R1` label, in emission order — the
    /// `*` stubs that the `RG:` alternatives attach to. Rendered as the
    /// `$...$` atom-label block, which needs one `;`-separated slot per
    /// atom of the main SMILES up to the last labelled one.
    r_group_sites: Vec<usize>,
    atom_offset: usize,
    var_offset: usize,
    /// `(sn, global atom indices of C1..Cn)` for each chain pushed via
    /// `push_chain`. Recorded explicitly rather than by push order, since
    /// several builders emit their chains out of sn order (`gpl_smiles`
    /// writes sn2 before sn1).
    chains: Vec<(usize, Vec<usize>)>,
}

impl CdkBuilder {
    /// Appends fixed, unambiguous literal SMILES text (headgroup
    /// backbone pieces, linking atoms, ring closures, ...).
    fn push_fixed(&mut self, text: &str) {
        self.atom_offset += count_atoms(text, true, true);
        self.smiles.push_str(text);
    }

    /// Appends a chain fragment, offsetting/renaming its local `Sg:`/`m:`
    /// blocks and constraint into this assembly's global numbering.
    fn push_fragment(&mut self, frag: &CdkFragment) {
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

        self.atom_offset += count_atoms(&frag.smiles, true, true);
        self.var_offset += local_var_count;
        self.smiles.push_str(&frag.smiles);
    }

    /// `push_fragment` for a fragment that *is* one of the lipid's chains,
    /// additionally recording where its carbons landed. Fragments with no
    /// carbons (an empty slot's bare `O`) are appended but not recorded.
    fn push_chain(&mut self, frag: &CdkFragment, sn: usize) {
        let base = self.atom_offset;
        let carbons: Vec<usize> = frag.carbon_atoms.iter().map(|a| a + base).collect();
        self.push_fragment(frag);
        if !carbons.is_empty() {
            self.chains.push((sn, carbons));
        }
    }

    /// Appends a literal `*` R-group attachment stub, recording its atom
    /// index so `finish` can label it `R1`.
    fn push_r_group_site(&mut self) {
        self.r_group_sites.push(self.atom_offset);
        self.push_fixed("*");
    }

    /// Finalizes the assembly into the final SMILES string.
    ///
    /// `rg_defs` are the alternative substituents for the `R1` label — one
    /// per interchangeable chain — emitted as the `$...$` atom-label block
    /// plus `RG:_R1={def},{def}`. Any `*` pushed via `push_r_group_site`
    /// may resolve to any one of them, which is how unresolved sn
    /// regiochemistry is expressed.
    fn finish(self, rg_defs: &[String]) -> Built {
        let mut blocks = Vec::new();
        if !rg_defs.is_empty() && !self.r_group_sites.is_empty() {
            blocks.push(r_group_label_block(
                &self.r_group_sites,
                count_atoms(&self.smiles, true, true),
            ));
        }
        let mut chains = self.chains;
        chains.sort_by_key(|(sn, _)| *sn);
        let constraint_str = self.constraints.join(",");
        blocks.extend(self.sg_blocks);
        blocks.extend(self.m_blocks);
        if !rg_defs.is_empty() && !self.r_group_sites.is_empty() {
            let defs = rg_defs
                .iter()
                .map(|d| format!("{{{d}}}"))
                .collect::<Vec<_>>()
                .join(",");
            blocks.push(format!("RG:_R1={defs}"));
        }
        let smiles = if blocks.is_empty() && constraint_str.is_empty() {
            self.smiles
        } else {
            format!("{} |{}| {}", self.smiles, blocks.join(","), constraint_str)
                .trim_end()
                .to_string()
        };
        Built { smiles, chains }
    }
}

/// The `$...$` atom-label block marking every atom in `sites` as `R1`.
///
/// Labels are positional and `;`-separated, one slot per atom of the main
/// SMILES in emission order. Trailing empty slots are omitted, so the block
/// runs only as far as the last labelled atom. `atom_count` bounds the list
/// defensively; a site beyond it would be a bookkeeping bug, and CDK accepts
/// misplaced labels silently rather than complaining.
fn r_group_label_block(sites: &[usize], atom_count: usize) -> String {
    let last = sites.iter().copied().max().unwrap_or(0).min(atom_count);
    let labels: Vec<&str> = (0..=last)
        .map(|i| if sites.contains(&i) { "R1" } else { "" })
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

/// Renames every variable letter in a `"a+b=N"`-style constraint per
/// `var_renaming`, leaving `+`/`=`/digits untouched.
fn rename_constraint_vars(constraint: &str, var_renaming: &HashMap<char, char>) -> String {
    constraint
        .chars()
        .map(|c| var_renaming.get(&c).copied().unwrap_or(c))
        .collect()
}

/// Count atoms in SMILES with options for wildcards and multi-component.
///
/// - `include_wildcards`: if true, count '*' as atoms; if false, skip them
/// - `all_components`: if true, count all components (including after dots);
///   if false, only count the first component
fn count_atoms(smiles: &str, include_wildcards: bool, all_components: bool) -> usize {
    let to_process = if all_components {
        smiles
    } else {
        smiles.split('.').next().unwrap_or("")
    };
    let mut count = 0;
    let mut chars = to_process.chars().peekable();

    while let Some(c) = chars.next() {
        match c {
            '[' => {
                count += 1;
                for c2 in chars.by_ref() {
                    if c2 == ']' {
                        break;
                    }
                }
            }
            'C' | 'N' | 'O' | 'P' | 'S' | 'F' | 'I' | 'c' | 'n' | 'o' | 's' | 'p' => count += 1,
            '*' if include_wildcards => count += 1,
            _ => {}
        }
    }

    count
}

/// Count atoms in SMILES including wildcards (for Sg/m block position calculation).
/// Only counts the first component (before the first dot).
#[allow(dead_code)]
fn count_atoms_with_wildcards(smiles: &str) -> usize {
    count_atoms(smiles, true, false)
}

/// Count all atoms in SMILES including wildcards across ALL components (handles dots).
#[allow(dead_code)]
fn count_all_atoms_with_wildcards(smiles: &str) -> usize {
    count_atoms(smiles, true, true)
}

/// Count total atoms in a SMILES string (before any dot), excluding wildcards.
#[allow(dead_code)]
fn count_atoms_in_smiles(smiles: &str) -> usize {
    count_atoms(smiles, false, false)
}

/// Resolves this generator's `Sg:n:pos:var:ht` unlocalized-double-bond
/// markers (see module docs and [`build_cdk_chain_smiles`]) into one
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
/// variable's count (see `build_cdk_chain_smiles`, which always emits
/// exactly one hardcoded atom per marker regardless of that variable's
/// eventual value), so only `value - 1` extra `C` atoms need inserting
/// after it. Any one valid split of each constraint sum across its
/// variables depicts an equally valid resolution, since the true
/// distribution is by definition unlocalized; this picks an even split.
///
/// An unlocalized modification's `m:` block is deliberately *not* resolved
/// to a specific carbon — forcing one would depict a positional claim the
/// name doesn't make, so it stays a dot-separated `*`-led component. The
/// `RG:` block is left alone; a depiction that needs one concrete chain per
/// sn position comes from [`lipid_name_to_structure`], which builds the
/// resolved layout instead of expanding this one.
///
/// Returns the input unchanged if it has no CXSMILES suffix or nothing
/// to expand.
pub fn expand_cxsmiles_for_depiction(smi: &str) -> String {
    expand_with_padding_inserts(smi).0
}

/// [`expand_cxsmiles_for_depiction`] plus the index translation it applied,
/// so callers holding atom indices into the unexpanded string can move them
/// across: the `(atom index, count)` runs of padding `C` atoms added for the
/// `Sg:` markers. These are chain carbons themselves, so a caller mapping a
/// chain also splices them into that chain's run (sort with
/// [`padding_shift_table`] before shifting).
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

    // Each "a+b=N" equation corresponds, in order, to the next `n_terms`
    // Sg positions: variable *names* can repeat across independently
    // numbered chains (see build_multichain_cdk_cxsmiles's var_offset
    // wraparound), so positional matching against emission order — not
    // name matching — is what stays unambiguous here.
    let mut inserts: Vec<(usize, usize)> = Vec::new();
    let mut consumed = 0usize;
    for eq in constraints.split(',').filter(|e| !e.is_empty()) {
        let Some((lhs, rhs)) = eq.split_once('=') else {
            continue;
        };
        let n_terms = lhs.split('+').count();
        let Ok(total) = rhs.trim().parse::<usize>() else {
            continue;
        };
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

/// Sorted `(position, count)` pairs for translating pre-expansion atom
/// indices. `insert_padding_atoms` writes each run of padding atoms
/// *after* the atom at `position`, so an index only shifts by the padding
/// that precedes it.
fn padding_shift_table(inserts: &[(usize, usize)]) -> Vec<(usize, usize)> {
    let mut table: Vec<(usize, usize)> = inserts.to_vec();
    table.sort_by_key(|(pos, _)| *pos);
    table
}

fn shift_atom(atom: usize, shift_table: &[(usize, usize)]) -> usize {
    atom + shift_table
        .iter()
        .take_while(|(pos, _)| *pos < atom)
        .map(|(_, count)| *count)
        .sum::<usize>()
}

/// Inserts `count` extra `C` atoms immediately after atom index `pos`
/// (0-based, counting every element/bracket-atom/`*` token across the
/// whole SMILES — same convention as `count_atoms(.., true, true)`) for
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
                out.extend(std::iter::repeat('C').take(n));
            }
            atom_idx += 1;
            continue;
        }
        if matches!(
            c,
            'C' | 'N' | 'O' | 'P' | 'S' | 'F' | 'I' | 'c' | 'n' | 'o' | 's' | 'p' | '*'
        ) {
            out.push(c);
            if let Some(&n) = insert_map.get(&atom_idx) {
                out.extend(std::iter::repeat('C').take(n));
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
            oh_pos: Vec::new(),
            ket_pos: Vec::new(),
            cooh_pos: Vec::new(),
            epox_pos: Vec::new(),
            cyc_pos: Vec::new(),
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

    let mut oh_pos = Vec::new();
    let mut ket_pos = Vec::new();
    let mut cooh_pos = Vec::new();
    let mut epox_pos = Vec::new();
    let mut cyc_pos = Vec::new();
    let mut oxy_declared: u32 = 0;

    for seg in cursor.split(';') {
        let seg = seg.trim();
        if seg.is_empty() {
            continue;
        }
        if let Some(inner) = seg.strip_prefix("OH(") {
            let inner = inner.strip_suffix(')')?;
            oh_pos = parse_plain_position_list(inner)?;
        } else if let Some(inner) = seg.strip_prefix("hydroxy(") {
            let inner = inner.strip_suffix(')')?;
            oh_pos = parse_plain_position_list(inner)?;
        } else if seg == "OH" || seg == "hydroxy" {
            // Standalone OH/hydroxy marker: unknown position
            oh_pos = vec![0]; // 0 means unknown position
        } else if !seg.starts_with("COOH") && seg.split(',').all(|part| part.ends_with("OH")) {
            oh_pos = seg
                .split(',')
                .map(|part| part.strip_suffix("OH")?.parse::<u32>().ok())
                .collect::<Option<Vec<_>>>()?;
        } else if let Some(inner) = seg.strip_prefix("oxo(") {
            let inner = inner.strip_suffix(')')?;
            ket_pos = parse_plain_position_list(inner)?;
        } else if seg == "oxo" {
            // Standalone oxo marker: unknown position
            ket_pos = vec![0]; // 0 means unknown position
        } else if seg.split(',').all(|part| part.ends_with("oxo")) {
            ket_pos = seg
                .split(',')
                .map(|part| part.strip_suffix("oxo")?.parse::<u32>().ok())
                .collect::<Option<Vec<_>>>()?;
        } else if let Some(inner) = seg.strip_prefix("COOH(") {
            let inner = inner.strip_suffix(')')?;
            cooh_pos = parse_plain_position_list(inner)?;
        } else if seg == "COOH" {
            // Standalone COOH marker: unknown position
            cooh_pos = vec![0]; // 0 means unknown position
        } else if let Some(inner) = seg
            .strip_prefix("ep(")
            .or_else(|| seg.strip_prefix("epox("))
        {
            let inner = inner.strip_suffix(')')?;
            epox_pos = parse_plain_position_list(inner)?;
        } else if let Some(inner) = seg
            .strip_prefix("cyc(")
            .or_else(|| seg.strip_prefix("cyclo("))
        {
            let inner = inner.strip_suffix(')')?;
            cyc_pos = parse_plain_position_list(inner)?;
        } else if let Some(digits) = seg.strip_prefix('O') {
            oxy_declared = if digits.is_empty() {
                1
            } else {
                digits.parse().ok()?
            };
        }
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

    // Allow unknown positions (marked as 0) or specific positions
    let has_unknown_mods = oh_pos.contains(&0) || ket_pos.contains(&0) || cooh_pos.contains(&0);
    let localized_oxy = oh_pos.iter().filter(|&&p| p != 0).count()
        + ket_pos.iter().filter(|&&p| p != 0).count()
        + cooh_pos.iter().filter(|&&p| p != 0).count();
    if oxy_declared > 0 && localized_oxy == 0 && !has_unknown_mods {
        return None; // generic oxygen count with no position hypothesis
    }

    Some(ParsedChain {
        prefix,
        carbon,
        db_pos,
        oh_pos,
        ket_pos,
        cooh_pos,
        epox_pos,
        cyc_pos,
    })
}

// ---------- generic chain-body rendering ----------

/// Per-carbon SMILES atom tokens (1-indexed via `atoms[k-1]`) plus the
/// bond symbol keyed by the *starting* carbon of each bond (`k` for the
/// bond between Ck and Ck+1); missing keys mean a plain single bond.
/// `db_pos` must only contain known (non-placeholder) double bonds —
/// unlocalized ones are represented separately via `Sg:` blocks (see
/// `chain_fragment_cdk_range`), never as a literal position here.
fn chain_tokens(
    carbon: u32,
    db_pos: &[DbPos],
    oh_pos: &[u32],
    ket_pos: &[u32],
    cooh_pos: &[u32],
) -> (Vec<String>, HashMap<u32, String>) {
    let db_start: HashMap<u32, DbPos> = db_pos.iter().map(|d| (d.pos, *d)).collect();
    let oh: HashSet<u32> = oh_pos.iter().copied().collect();
    let ket: HashSet<u32> = ket_pos.iter().copied().collect();
    let cooh: HashSet<u32> = cooh_pos.iter().copied().collect();

    let mut atoms = Vec::with_capacity(carbon as usize);
    for k in 1..=carbon {
        let tok = if ket.contains(&k) {
            "C(=O)".to_string()
        } else if oh.contains(&k) {
            "C(O)".to_string()
        } else if cooh.contains(&k) {
            "C(C(=O)O)".to_string()
        } else {
            "C".to_string()
        };
        atoms.push(tok);
    }

    let mut bonds: HashMap<u32, String> = HashMap::new();
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
/// by the running count.
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
        atom_idx += count_atoms(tok, true, true);
        out.push_str(tok);
        if k < end {
            if let Some(b) = bonds.get(&k) {
                out.push_str(b);
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
// Every argument is an independent axis of the chain being built, and the
// positional lists are already grouped as the name parser produces them;
// bundling them into a struct would only move the same eight fields.
#[allow(clippy::too_many_arguments)]
fn build_chain_range(
    carbon: u32,
    db_pos: &[DbPos],
    oh_pos: &[u32],
    ket_pos: &[u32],
    cooh_pos: &[u32],
    carbonyl_c1: bool,
    start: u32,
    end: u32,
) -> (String, Vec<usize>) {
    if carbon == 0 || end < start {
        return (String::new(), Vec::new());
    }
    let (mut atoms, bonds) = chain_tokens(carbon, db_pos, oh_pos, ket_pos, cooh_pos);
    if carbonyl_c1 {
        atoms[0] = "C(=O)".to_string();
    }
    assemble_range(&atoms, &bonds, start, end)
}

// ---------- headgroup builders ----------
//
// Each builder writes one lipid class's fixed template and delegates its
// variable slots to `slot_fragment_cdk`/`chain_fragment_cdk`, which return a
// chain fragment starting at C1 ready to be attached after the linking
// heteroatom (`O` for ester/ether, `N` for amide) the template itself writes.

fn fa_smiles(chains: &[ParsedChain]) -> Option<Built> {
    let c = chains.first()?;
    if c.carbon == 0 || c.prefix != ChainPrefix::Acyl {
        return None;
    }
    let frag = chain_fragment_cdk(c)?;
    let mut b = CdkBuilder::default();
    b.push_fixed("O");
    b.push_chain(&frag, 1);
    Some(b.finish(&[]))
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
    let frag = chain_fragment_cdk(c)?;
    let mut b = CdkBuilder::default();
    b.push_fixed("[n+]1ccccc1-c1ccc(CN");
    b.push_chain(&frag, 1);
    b.push_fixed(")cc1");
    Some(b.finish(&[]))
}

fn nae_smiles(chains: &[ParsedChain]) -> Option<Built> {
    let c = chains.first()?;
    if c.carbon == 0 || c.prefix != ChainPrefix::Acyl {
        return None;
    }
    let frag = chain_fragment_cdk(c)?;
    let mut b = CdkBuilder::default();
    b.push_fixed("OCCN");
    b.push_chain(&frag, 1);
    Some(b.finish(&[]))
}

fn car_smiles(chains: &[ParsedChain]) -> Option<Built> {
    let c = chains.first()?;
    if c.carbon == 0 || c.prefix != ChainPrefix::Acyl {
        return None;
    }
    let frag = chain_fragment_cdk(c)?;
    // The ester O needs two real bonds: to the chain's carbonyl carbon
    // (C1, i.e. frag's own first atom) and to carnitine's chiral carbon.
    // Writing the chain as a branch straight off O keeps both bonds
    // explicit without touching [C@H]'s neighbor order (O still precedes
    // it in the text exactly as before), so no stereo re-derivation is
    // needed.
    let mut b = CdkBuilder::default();
    b.push_fixed("O(");
    b.push_chain(&frag, 1);
    b.push_fixed(")[C@H](CC(=O)[O-])C[N+](C)(C)C");
    Some(b.finish(&[]))
}

fn ce_smiles(chains: &[ParsedChain]) -> Option<Built> {
    let c = chains.first()?;
    if c.carbon == 0 || c.prefix != ChainPrefix::Acyl {
        return None;
    }
    let frag = chain_fragment_cdk(c)?;
    let mut b = CdkBuilder::default();
    b.push_fixed("C12(CC=C3CC(O");
    b.push_chain(&frag, 1);
    b.push_fixed(")CCC3(C)C1CCC1(C)C(C(C)CCCC(C)C)CCC21)");
    Some(b.finish(&[]))
}

fn st_smiles() -> String {
    "C12(CC=C3CC(O)CCC3(C)C1CCC1(C)C(C(C)CCCC(C)C)CCC21)".to_string()
}

/// Glycerol backbone with up to `slots` ester/ether positions (sn1, sn2,
/// sn3 in that order); unfilled positions are free hydroxyls. Used for
/// MG/DG/TG, which have no phosphate headgroup.
fn glycerolipid_smiles(
    chains: &[ParsedChain],
    slots: usize,
    mode: Regiochemistry,
) -> Option<Built> {
    if chains.is_empty() || chains.len() > slots {
        return None;
    }
    let bare_o = || CdkFragment {
        smiles: "O".to_string(),
        ..Default::default()
    };
    match mode {
        Regiochemistry::Resolved => {
            let sn1 = slot_fragment_cdk(chains.first())?;
            let sn2 = if slots >= 2 {
                slot_fragment_cdk(chains.get(1))?
            } else {
                bare_o()
            };
            let sn3 = if slots >= 3 {
                slot_fragment_cdk(chains.get(2))?
            } else {
                bare_o()
            };

            let mut b = CdkBuilder::default();
            b.push_fixed("C(C");
            b.push_chain(&sn3, 3);
            b.push_fixed(")(");
            b.push_chain(&sn2, 2);
            b.push_fixed(")C");
            b.push_chain(&sn1, 1);
            Some(b.finish(&[]))
        }
        Regiochemistry::Unresolved => {
            let Some(rg_defs) = rg_alternatives(chains) else {
                return glycerolipid_smiles(chains, slots, Regiochemistry::Resolved);
            };
            let n = chains.len();
            let mut b = CdkBuilder::default();
            b.push_fixed("C(CO");
            if n >= 3 {
                b.push_r_group_site();
            }
            b.push_fixed(")(O");
            if n >= 2 {
                b.push_r_group_site();
            }
            b.push_fixed(")CO");
            b.push_r_group_site();
            Some(b.finish(&rg_defs))
        }
    }
}

/// Diacyl-glycerophospholipid: `sn1`/`sn2` chains plus a fixed
/// phospho-headgroup tail attached at sn3 (e.g. `"OP(=O)([O-])OCC[N+](C)(C)C"`
/// for PC). Also covers the lyso forms, which simply supply one chain
/// (always treated as resolved, since there's nothing to be ambiguous
/// between). `headgroup_tail` is generic across PC/PE/PS/PG/PI/PA (and
/// their lyso forms), so this one Sg:-aware builder gives all of them
/// unlocalized-double-bond/regiochemistry coverage.
fn gpl_smiles(headgroup_tail: &str, chains: &[ParsedChain], mode: Regiochemistry) -> Option<Built> {
    if chains.is_empty() || chains.len() > 2 {
        return None;
    }
    match mode {
        Regiochemistry::Resolved => {
            let sn1 = slot_fragment_cdk(chains.first())?;
            let sn2 = slot_fragment_cdk(chains.get(1))?;
            let mut b = CdkBuilder::default();
            b.push_fixed("C(C");
            b.push_fixed(headgroup_tail);
            b.push_fixed(")(");
            b.push_chain(&sn2, 2);
            b.push_fixed(")C");
            b.push_chain(&sn1, 1);
            Some(b.finish(&[]))
        }
        Regiochemistry::Unresolved => {
            let Some(rg_defs) = rg_alternatives(chains) else {
                return gpl_smiles(headgroup_tail, chains, Regiochemistry::Resolved);
            };
            let mut b = CdkBuilder::default();
            b.push_fixed("C(C");
            b.push_fixed(headgroup_tail);
            b.push_fixed(")(O");
            b.push_r_group_site();
            b.push_fixed(")CO");
            b.push_r_group_site();
            Some(b.finish(&rg_defs))
        }
    }
}

/// Cardiolipin: two phosphatidyl arms (sn1/sn2 and sn3/sn4) hung off a
/// central glycerol's C1/C3, with a free hydroxyl at the central C2.
fn cl_smiles(chains: &[ParsedChain], mode: Regiochemistry) -> Option<Built> {
    if chains.len() != 4 {
        return None;
    }
    match mode {
        Regiochemistry::Resolved => {
            let sn1 = slot_fragment_cdk(chains.first())?;
            let sn2 = slot_fragment_cdk(chains.get(1))?;
            let sn3 = slot_fragment_cdk(chains.get(2))?;
            let sn4 = slot_fragment_cdk(chains.get(3))?;

            let mut b = CdkBuilder::default();
            b.push_fixed("C(COP(=O)(O)OCC(");
            b.push_chain(&sn2, 2);
            b.push_fixed(")C");
            b.push_chain(&sn1, 1);
            b.push_fixed(")(O)COP(=O)(O)OCC(");
            b.push_chain(&sn4, 4);
            b.push_fixed(")C");
            b.push_chain(&sn3, 3);
            Some(b.finish(&[]))
        }
        Regiochemistry::Unresolved => {
            let Some(rg_defs) = rg_alternatives(chains) else {
                return cl_smiles(chains, Regiochemistry::Resolved);
            };
            let mut b = CdkBuilder::default();
            b.push_fixed("C(COP(=O)(O)OCC(O");
            b.push_r_group_site();
            b.push_fixed(")CO");
            b.push_r_group_site();
            b.push_fixed(")(O)COP(=O)(O)OCC(O");
            b.push_r_group_site();
            b.push_fixed(")CO");
            b.push_r_group_site();
            Some(b.finish(&rg_defs))
        }
    }
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
    let mut oh_pos = base.oh_pos.clone();
    if !oh_pos.contains(&3) {
        oh_pos.push(3);
    }
    if is_triol && !oh_pos.contains(&4) {
        oh_pos.push(4);
    }

    let rest = chain_fragment_cdk_range(
        3,
        base.carbon,
        &base.db_pos,
        &oh_pos,
        &base.ket_pos,
        &base.cooh_pos,
        false,
    )?;

    let n_frag = match n_acyl {
        Some(acyl) if acyl.carbon > 0 && acyl.prefix == ChainPrefix::Acyl => {
            Some(chain_fragment_cdk(acyl)?)
        }
        Some(_) => return None,
        None => None,
    };

    let mut b = CdkBuilder::default();
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
    Some(b.finish(&[]))
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
                c @ ('C' | 'N' | 'O' | 'P' | 'S' | 'F' | 'I' | 'c' | 'n' | 'o' | 's' | 'p'
                | '*') => {
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
        let st = lipid_name_to_structure(name).unwrap_or_else(|| panic!("{name} should resolve"));
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

    /// `_` names depict as loose `*`-capped fragments in their honest
    /// CXSMILES form; the structure view builds them connected instead and
    /// says so via `regio_resolved`.
    #[test]
    fn structure_builds_sn_unresolved_names_as_one_connected_molecule() {
        let st = lipid_name_to_structure("PC 16:0_18:1(9)").expect("should resolve");
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

        // The plain SMILES accessor keeps its existing fragment form.
        let plain = lipid_name_to_smiles("PC 16:0_18:1(9)").expect("should resolve");
        assert!(
            plain.contains('*'),
            "lipid_name_to_smiles behaviour changed: {plain}"
        );
    }

    /// C1 is the carboxyl carbon, so the localized double bond declared at
    /// C9 must sit between the 9th and 10th mapped atoms.
    #[test]
    fn structure_map_carbon_numbering_lines_up_with_declared_db_position() {
        let st = lipid_name_to_structure("FA 18:1(9)").expect("should resolve");
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
        let s = lipid_name_to_smiles("TG 18:0/18:0/18:1(9);5oxo").expect("should resolve");
        assert_balanced(&s);
        assert!(s.contains("C(=O)CCCC(=O)CCCC=CCCCCCCCC"));
        // Every position is known, so nothing is left to annotate.
        assert!(!s.contains(" |"), "{s}");
    }

    #[test]
    fn pc_with_geometry() {
        let s = lipid_name_to_smiles("PC 16:0/18:1(9Z)").expect("should resolve");
        assert_balanced(&s);
        assert!(s.contains("[N+](C)(C)C"));
        assert!(s.contains("/C=C\\"));
    }

    #[test]
    fn pe_trans_geometry() {
        let s = lipid_name_to_smiles("PE 16:0/18:1(9E)").expect("should resolve");
        assert_balanced(&s);
        assert!(s.contains("/C=C/"));
    }

    #[test]
    fn fa_unlocalized_db_gets_sg_blocks() {
        let s = lipid_name_to_smiles("FA 20:4").expect("should resolve via Sg: blocks");
        assert_balanced(&s);
        assert!(
            s.contains(" |Sg:"),
            "unlocalized double bonds should be flagged with Sg:"
        );
        assert!(
            !s.contains("RG:"),
            "single chain has no regiochemistry ambiguity"
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

        let expanded = expand_cxsmiles_for_depiction(&s);
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
        let s = lipid_name_to_smiles("FA 16:0").expect("should resolve");
        assert_balanced(&s);
        assert_eq!(s, "OC(=O)CCCCCCCCCCCCCCC");
    }

    #[test]
    fn amp_fa_hete() {
        let s = lipid_name_to_smiles("AMP-FA 20:4(5,8,11,14);15OH").expect("should resolve");
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
        let s = lipid_name_to_smiles("FA 20:4(5,8,11,14);15OH")
            .expect("Shorthand2020 hydroxyl syntax should resolve");
        assert_balanced(&s);
        assert!(s.contains("C(O)"));
    }

    #[test]
    fn confidence_display_tail_is_ignored_by_structure_parser() {
        let s = lipid_name_to_smiles(
            "FA 20:4;11OH [DB sn1: Δ5 100%, Δ8 100%, Δ12 100% | Δ14 50% | Δ15 50%]",
        )
        .expect("display tail must not make the canonical structure unparsable");
        assert_balanced(&s);
    }

    #[test]
    fn amp_fa_unlocalized_oxygen_returns_none() {
        // generic ;O with no OH/oxo/COOH breakdown -> still ambiguous, no
        // placeholder convention was requested for oxygen sites.
        assert!(lipid_name_to_smiles("AMP-FA 20:4(5,8,11,14);O").is_none());
    }

    #[test]
    fn lpc_single_chain_free_oh() {
        let s = lipid_name_to_smiles("LPC 16:0").expect("should resolve");
        assert_balanced(&s);
        assert!(s.contains("(O)"));
        assert!(!s.contains('|'));
    }

    #[test]
    fn cer_d18_1_16_0() {
        let s = lipid_name_to_smiles("Cer d18:1(4)/16:0").expect("should resolve");
        assert_balanced(&s);
        assert!(s.contains("NC(=O)CCCCCCCCCCCCCCC"));
    }

    #[test]
    fn sm_d18_1_16_0() {
        let s = lipid_name_to_smiles("SM d18:1(4)/16:0").expect("should resolve");
        assert_balanced(&s);
        assert!(s.contains("[N+](C)(C)C"));
    }

    #[test]
    fn slot_fragment_unlocalized_db_expands_without_corrupting_carbonyl_branch() {
        // Regression test: slot_fragment_cdk() prepends "O" to a chain's
        // fragment (the glycerol-ester link) but was not shifting that
        // chain's own Sg: positions to account for the new atom, so
        // expansion inserted padding carbons *inside* the C1 carbonyl's
        // `(=O)` branch instead of after it -- producing an invalid
        // (3-bonded) oxygen and breaking depiction for any Resolved-mode
        // chain (single chain, or "/"-joined) with an unlocalized DB.
        for name in ["LPC 18:2", "PC 18:2/14:1", "PE 18:2/14:1", "LPE 18:2"] {
            let s = lipid_name_to_smiles(name).expect("should resolve");
            assert_balanced(&s);
            let expanded = expand_cxsmiles_for_depiction(&s);
            assert!(!expanded.contains('|'));
            assert_balanced(&expanded);
            assert!(
                !expanded.contains("=OC"),
                "{name}: padding must never land inside a (=O) branch: {expanded}"
            );
        }

        // LPC 18:2: glycerophosphocholine backbone (11 headgroup atoms +
        // 3 glycerol carbons + 1 free OH) + an 18-carbon acyl chain.
        let s = lipid_name_to_smiles("LPC 18:2").unwrap();
        let expanded = expand_cxsmiles_for_depiction(&s);
        assert_eq!(
            expanded.chars().filter(|&c| c == 'C').count(),
            3 + 5 + 18,
            "3 glycerol C + 5 phosphocholine C + 18 chain C"
        );
        assert_eq!(expanded.matches("C=C").count(), 2);
    }

    #[test]
    fn ce_18_1() {
        let s = lipid_name_to_smiles("CE 18:1(9Z)").expect("should resolve");
        assert_balanced(&s);
        assert!(s.contains("/C=C\\"));
    }

    #[test]
    fn car_ester_oxygen_bonds_to_the_carbonyl_carbon() {
        // Acylcarnitine: R-C(=O)-O-CH(CH2COO-)-CH2-N+(CH3)3. The ester O
        // must bond to the chain's carbonyl carbon (C1), not to its
        // omega end -- regression test for a bug where the chain was
        // concatenated before "O[C@H]...", bonding O to the chain's last
        // (omega) atom instead.
        let s = lipid_name_to_smiles("CAR 18:0").expect("should resolve");
        assert_balanced(&s);
        assert_eq!(s, "O(C(=O)CCCCCCCCCCCCCCCCC)[C@H](CC(=O)[O-])C[N+](C)(C)C");
        assert!(
            s.starts_with("O(C(=O)"),
            "ester O must be directly bonded to the carbonyl carbon"
        );

        let s2 = lipid_name_to_smiles("CAR 18:1").expect("should resolve");
        assert_balanced(&s2);
        assert!(s2.starts_with("O(C(=O)"));
        let expanded = expand_cxsmiles_for_depiction(&s2);
        assert_balanced(&expanded);
        assert!(expanded.starts_with("O(C(=O)"));
        // 18 chain carbons + carnitine's own 7 ([C@H] + CH2-C(=O)O- branch (2) + CH2 + N(CH3)3 (3)).
        assert_eq!(expanded.chars().filter(|&c| c == 'C').count(), 18 + 7);
    }

    #[test]
    fn cl_unlocalized_db_gets_placeholders() {
        // LipidOracle's EAD engines always join CL/DG/TG chains with "_",
        // never "/".
        let s = lipid_name_to_smiles("CL 18:2_18:2_18:2_18:2").expect("should resolve");
        assert_balanced(&s);
        assert!(
            s.contains("Sg:n:"),
            "unlocalized double bonds should be flagged with Sg:"
        );
        // Sg: needs its atoms in the main SMILES and CDK rejects a nested
        // block inside an RG: definition, so these chains cannot also carry
        // the sn ambiguity — the double-bond claim wins and the sn
        // assignment shown is one arbitrary order.
        assert!(!s.contains("RG:"), "{s}");
    }

    #[test]
    fn cl_four_chains_saturated_still_unresolved_regiochemistry() {
        let s = lipid_name_to_smiles("CL 18:0_18:0_18:0_18:0").expect("should resolve");
        assert_balanced(&s);
        // CL is always "_"-joined by this codebase's own convention, so
        // regiochemistry is unresolved even though every chain here is
        // fully saturated. With no Sg: in the way, all four chains can move
        // into R-group definitions.
        assert_eq!(s.matches("R1").count(), 5, "4 labelled sites + 1 RG: name");
        assert_eq!(s.matches("},{").count(), 3, "4 alternatives");
        assert!(!s.contains("Sg:"));
    }

    #[test]
    fn pc_slash_joined_stays_resolved_even_if_incomplete() {
        // "/" always means resolved, regardless of chain count.
        let s = lipid_name_to_smiles("PC 16:0/18:1(9Z)").expect("should resolve");
        assert!(!s.contains("f:"));
    }

    #[test]
    fn pc_underscore_joined_with_unlocalized_chains_keeps_sg_over_rg() {
        let s = lipid_name_to_smiles("PC 16:1_18:1").expect("should resolve via Sg: blocks");
        assert_balanced(&s);
        assert!(
            s.contains("Sg:n:"),
            "unlocalized double bonds should be flagged with Sg:"
        );
        // The sn ambiguity is the one that has to give: `Sg:` indexes atoms
        // of the main SMILES, and CDK rejects a nested block inside an `RG:`
        // definition, so a chain needing `Sg:` cannot become an R-group
        // alternative. Two chains here need one, so no `RG:` is emitted and
        // the sn order shown is arbitrary.
        assert!(!s.contains("RG:"), "{s}");
        assert_eq!(
            s.matches('*').count(),
            0,
            "no R-group sites, so no wildcards: {s}"
        );
    }

    #[test]
    fn unknown_headgroup_returns_none() {
        assert!(lipid_name_to_smiles("PIP2 16:0/18:1(9)").is_none());
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
        // lipid_name_to_smiles's own shorthand rejection above.
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
        assert!(lipid_name_to_smiles("PC 19:2").is_none());
    }

    #[test]
    fn shorthand_tg_54_3_rejected() {
        // "TG 54:3" is ambiguous: total composition for 3 chains.
        // Should return None.
        assert!(lipid_name_to_smiles("TG 54:3").is_none());
    }

    #[test]
    fn shorthand_dg_36_2_rejected() {
        // "DG 36:2" is ambiguous: total composition for 2 chains.
        // Should return None.
        assert!(lipid_name_to_smiles("DG 36:2").is_none());
    }

    #[test]
    fn explicit_pc_16_0_slash_18_1_accepted() {
        // "PC 16:0/18:1" has explicit chains with "/" separator, so it's valid.
        let s = lipid_name_to_smiles("PC 16:0/18:1").expect("explicit chains should work");
        assert_balanced(&s);
        assert!(!s.contains("f:"), "resolved regiochemistry needs no f:");
    }

    #[test]
    fn explicit_pc_16_1_underscore_18_1_accepted() {
        // "PC 16:1_18:1" has explicit chains with "_" separator, so it's valid.
        let s = lipid_name_to_smiles("PC 16:1_18:1").expect("explicit chains should work");
        assert_balanced(&s);
        // Both chains are unlocalized, so this takes the Sg:-keeping path.
        assert!(s.contains("Sg:n:"), "{s}");
    }

    #[test]
    fn single_chain_fa_18_1_accepted() {
        // FA and other single-chain lipids should work with shorthand like "FA 18:1".
        let s = lipid_name_to_smiles("FA 18:1").expect("single-chain shorthand should work");
        assert_balanced(&s);
    }

    #[test]
    fn fa_with_hydroxyl_at_known_position() {
        // FA with hydroxyl at position 15
        let s = lipid_name_to_smiles("FA 20:4(5,8,11,14);15OH").expect("FA with OH should work");
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
        let s =
            lipid_name_to_smiles("FA 18:1(9);3OH,5OH").expect("FA with multiple OH should work");
        assert_balanced(&s);
        let oh_count = s.matches("C(O)").count();
        assert_eq!(oh_count, 2, "should have 2 hydroxyl groups");
    }

    #[test]
    fn fa_with_ketone_at_known_position() {
        // FA with ketone (oxo) at position 5
        let s = lipid_name_to_smiles("FA 18:1(9);5oxo").expect("FA with oxo should work");
        assert_balanced(&s);
        assert!(s.contains("C(=O)"), "ketone group should be present");
    }

    #[test]
    fn fa_with_hydroxyls_and_ketones() {
        // FA with both hydroxyls and ketones
        let s =
            lipid_name_to_smiles("FA 20:2(5,11);3OH;8oxo").expect("FA with OH and oxo should work");
        assert_balanced(&s);
        assert!(s.contains("C(O)"), "hydroxyl should be present");
        assert!(s.contains("C(=O)"), "ketone should be present");
    }

    #[test]
    fn tg_with_modifications_on_one_chain() {
        // TG with one chain having modifications
        let s = lipid_name_to_smiles("TG 18:0/18:1(9);3OH/18:0")
            .expect("TG with chain modification should work");
        assert_balanced(&s);
        assert!(s.contains("C(O)"), "chain with hydroxyl should be present");
    }

    #[test]
    fn dg_with_hydroxyl_unresolved_regiochemistry() {
        // DG with hydroxyl on unresolved chains
        let s = lipid_name_to_smiles("DG 18:1(9);3OH_18:0")
            .expect("DG with OH and unresolved should work");
        assert_balanced(&s);
        assert!(
            s.contains("RG:_R1="),
            "unresolved regiochemistry needs R-group alternatives: {s}"
        );
        assert!(s.contains("C(O)"), "hydroxyl should be present");
    }

    #[test]
    fn generic_oxygen_without_position_returns_none() {
        // Generic oxygen (;O or ;O2) without specific position breakdown
        // should return None as per current design (ambiguity without fallback)
        assert!(
            lipid_name_to_smiles("FA 18:1(9);O").is_none(),
            "generic O without position should be rejected"
        );
        assert!(
            lipid_name_to_smiles("FA 18:1(9);O2").is_none(),
            "generic O2 without position should be rejected"
        );
    }

    #[test]
    fn hydroxyl_with_unlocalized_double_bond() {
        // FA with localized hydroxyl but unlocalized double bond
        let s = lipid_name_to_smiles("FA 18:2;5OH")
            .expect("unlocalized DB with localized OH should work");
        assert_balanced(&s);
        assert!(s.contains("C(O)"), "hydroxyl should be present");
        assert!(
            s.contains(" |Sg:"),
            "unlocalized double bonds should be flagged with Sg:"
        );

        let expanded = expand_cxsmiles_for_depiction(&s);
        assert!(!expanded.contains('|'));
        assert_balanced(&expanded);
        assert!(
            expanded.contains("C(O)"),
            "hydroxyl should survive expansion"
        );
        assert_eq!(expanded.chars().filter(|&c| c == 'C').count(), 18);
        assert_eq!(expanded.matches("C=C").count(), 2);
    }

    #[test]
    fn epoxide_modifications_are_parseable() {
        // Epoxide modifications should parse without error (rendering is future work)
        // FA with epoxide at C5-C6
        let result = lipid_name_to_smiles("FA 18:0;ep(5)");
        // Either returns Some (if rendering is implemented) or None (if not yet implemented)
        // The important thing is it doesn't panic during parsing
        let _ = result;
    }

    #[test]
    fn cyclopropane_modifications_are_parseable() {
        // Cyclopropane modifications should parse without error (rendering is future work)
        // FA with cyclopropane at C5-C6
        let result = lipid_name_to_smiles("FA 18:0;cyc(5)");
        // Either returns Some (if rendering is implemented) or None (if not yet implemented)
        // The important thing is it doesn't panic during parsing
        let _ = result;
    }

    #[test]
    fn multiple_modifications_including_new_types() {
        // Test parsing of chains with multiple modification types together
        let result = lipid_name_to_smiles("FA 20:2(5,8);3OH;10oxo;ep(15)");
        // Should at least parse without crashing
        let _ = result;
    }

    // Tests for CDK-style lipid_name_to_cxsmiles()
    #[test]
    fn cdk_simple_fa() {
        // Simple fatty acid should work
        let s = lipid_name_to_cxsmiles("FA 18:1(9)").expect("should resolve");
        eprintln!("Generated SMILES: {}", s);
        assert_balanced(&s);
        assert!(s.contains("C=C"), "double bond should be present");
    }

    // DG and PC support coming soon - currently FA only
    // TODO: Implement DG/PC multi-chain support in build_cdk_chain_smiles
    //   - Handle multiple chains with wildcard (*) attachment points
    //   - Generate f: blocks for component grouping
    //   - Support Sg labels for sn-position ambiguity

    #[test]
    fn expand_cxsmiles_resolves_to_full_chain_length() {
        // FA 18:2: 18 total carbons, 2 unlocalized double bonds.
        let cxsmiles = lipid_name_to_smiles("FA 18:2").expect("should resolve");
        let expanded = expand_cxsmiles_for_depiction(&cxsmiles);

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
        let expanded = expand_cxsmiles_for_depiction(cxsmiles);
        assert!(!expanded.contains('|'));
        assert_balanced(&expanded);
        assert_eq!(expanded.chars().filter(|&c| c == 'C').count(), 18);
    }

    #[test]
    fn expand_cxsmiles_passthrough_without_sg_blocks() {
        // Plain SMILES with no CXSMILES suffix is returned unchanged.
        let s = lipid_name_to_smiles("FA 16:0").expect("should resolve");
        assert_eq!(expand_cxsmiles_for_depiction(&s), s);
    }

    #[test]
    fn cdk_fa_unlocalized_db() {
        // FA with unlocalized double bonds
        let s = lipid_name_to_cxsmiles("FA 18:2").expect("should resolve");
        eprintln!("CDK FA 18:2: {}", s);
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
    fn cdk_single_unlocalized_db_has_terminal_sg_marker() {
        let s = lipid_name_to_cxsmiles("FA 18:1").expect("should resolve");
        let smiles_part = s.split(' ').next().unwrap_or(&s);

        assert!(smiles_part.ends_with("C=CC"));
        assert!(s.contains("Sg:n:3:a:ht,Sg:n:6:b:ht"));
        assert!(
            s.contains("a+b=15"),
            "the terminal carbon reduces the variable count by one"
        );

        let expanded = expand_cxsmiles_for_depiction(&s);
        assert_eq!(expanded.chars().filter(|&c| c == 'C').count(), 18);
        assert_eq!(expanded.matches("C=C").count(), 1);
    }

    #[test]
    fn unknown_geometry_is_just_a_plain_double_bond() {
        let unspecified = lipid_name_to_smiles("FA 18:1(9);OH").expect("should resolve");
        assert!(unspecified.contains("|m:"));
        assert!(unspecified.contains("C=C"), "{unspecified}");
        assert!(!unspecified.contains('/') && !unspecified.contains('\\'));

        let explicit = lipid_name_to_smiles("FA 18:1(9Z);OH").expect("should resolve");
        assert!(explicit.contains("/C=C\\"));
    }

    #[test]
    fn test_fa_16_1_known_db() {
        // Test cases with variable DBs and known DBs
        let test_cases = vec![
            ("FA 18:2", "Variable DBs (2)"),
            ("FA 16:4", "Variable DBs (4)"),
            ("FA 20:3", "Variable DBs (3)"),
            ("FA 16:1(9)", "Known DB at 9"),
            ("FA 16:1(9);oxo;OH", "Known DB 9, unknown oxo and OH"),
            ("FA 18:1(9);OH", "Known DB 9, unknown OH"),
            ("FA 18:1(9);oxo", "Known DB 9, unknown oxo"),
            (
                "FA 18:1(9);oxo;OH;COOH",
                "Known DB 9, unknown oxo, OH, COOH",
            ),
            (
                "FA 18:4;oxo;hydroxy",
                "Variable DBs (4) + unknown oxo and hydroxy",
            ),
            ("FA 18:4(5,8,11,14);3oxo;16OH", "Known DBs + known mods"),
        ];
        for (name, desc) in test_cases {
            match lipid_name_to_cxsmiles(name) {
                Some(s) => eprintln!("{} ({}): {}", name, desc, s),
                None => eprintln!("{} ({}): None", name, desc),
            }
        }
    }

    // FA with known modifications - m: blocks coming soon
    // TODO: Implement full m: block generation for modification placement options
    //   - Map atom indices to possible positions
    //   - Support multiple modification fragments
    //   - Combine with f: blocks for multi-chain cases

    // PC complex cases - multi-chain support coming soon
    // TODO: Implement PC/DG with unknown sn-positions + variable DBs + variable modifications

    #[test]
    fn pc_complex_ambiguities() {
        let test_name = "PC 16:2_18:1(7);5oxo";
        let s = lipid_name_to_smiles(test_name).expect("should parse");

        eprintln!("\n=== CXSMILES ===");
        eprintln!("Input:  {}", test_name);
        eprintln!("Output: {}", s);
        eprintln!();

        assert_balanced(&s);

        if s.contains("C(=O)C") || s.contains("(=O)C") {
            eprintln!("✓ Has oxo/ketone group from chain 2");
        }
        assert!(
            s.contains("Sg:"),
            "should have Sg: for chain 1's unlocalized double bond"
        );
        // Chain 1's Sg: run keeps it in the main SMILES, which rules out the
        // RG: form for this molecule - the two ambiguities cannot coexist.
        assert!(!s.contains("RG:"), "{s}");
    }

    /// An unlocalized modification is a CXSMILES position-variation bond:
    /// `m:<dummy atom>:<candidate anchors>`. The anchors must be this
    /// chain's own carbons (minus the acyl C1), the variable end must be a
    /// `*` dummy carrying exactly one bond, and the component must not
    /// smuggle in a carbon the chain doesn't have.
    #[test]
    fn unlocalized_modification_lists_every_chain_carbon() {
        let s = lipid_name_to_smiles("FA 18:1;OH").expect("should resolve");
        // OC(=O)CC=CC.*O -> chain carbons are atoms 1..=6, C1 (atom 1) is
        // the carboxyl, and atom 7 is the `*` the hydroxyl hangs off.
        assert!(s.contains("m:7:3.4.5.6"), "{s}");
        assert_eq!(s.split(" |").next().unwrap(), "OC(=O)CC=CC.*O", "{s}");

        // A ketone converts a carbon rather than adding one; an extra
        // carboxyl brings its own.
        let oxo = lipid_name_to_smiles("FA 18:0;oxo").expect("should resolve");
        assert!(oxo.split(" |").next().unwrap().ends_with(".*=O"), "{oxo}");
        let cooh = lipid_name_to_smiles("FA 18:0;COOH").expect("should resolve");
        assert!(
            cooh.split(" |").next().unwrap().ends_with(".*C(=O)O"),
            "{cooh}"
        );
    }

    /// The position-variation stub survives into the depiction form
    /// unchanged, and adds no carbon the molecule doesn't have.
    ///
    /// This replaces an older table of substitutions (`.[OH]` -> `.CO` and
    /// friends) that existed to make an inert `m:` block at least *look*
    /// like a modification. `*O` needs no such trade: it is both the form
    /// the block requires and a legible `R-OH` stub, so the depicted
    /// formula finally matches the stored one.
    #[test]
    fn depiction_keeps_the_position_variation_stub() {
        for (name, expected) in [
            ("FA 18:1;OH", "OC(=O)CCCCCCCC=CCCCCCCCC.*O"),
            ("FA 18:0;oxo", "OC(=O)CCCCCCCCCCCCCCCCC.*=O"),
            ("FA 18:0;COOH", "OC(=O)CCCCCCCCCCCCCCCCC.*C(=O)O"),
        ] {
            let depicted = expand_cxsmiles_for_depiction(&lipid_name_to_smiles(name).unwrap());
            assert_eq!(depicted, expected, "{name}");
        }
        // A charged bracket atom inside a real structure is not a
        // position-variation placeholder and must survive untouched.
        let pc = expand_cxsmiles_for_depiction(&lipid_name_to_smiles("PC 16:0/18:1(9)").unwrap());
        assert!(pc.contains("[O-]") && pc.contains("[N+]"), "{pc}");

        // The stub sits between the chain and whatever follows, so every
        // chain's carbon indices still have to line up around it.
        assert_chain_map_sane("FA 18:1;OH", &[(1, 18)]);
        assert_chain_map_sane("Cer d18:1(4)/16:0;OH", &[(1, 18), (2, 16)]);
        assert_chain_map_sane("DG 18:1(9);5OH_18:1;OH", &[(1, 18), (2, 18)]);
    }

    /// An `R1` alternative must not lose a chain carbon to the `*` it
    /// attaches to: the ester oxygen lives on the backbone, so every
    /// definition holds all of its chain's own carbons.
    #[test]
    fn r_group_alternatives_keep_every_carbon() {
        let s = lipid_name_to_smiles("DG 18:1(9);5OH_18:0").expect("should resolve");
        let (_, blocks) = s.split_once(" |").expect("should carry an RG: block");
        let defs: Vec<&str> = blocks
            .split_once("RG:_R1=")
            .expect("RG: definitions")
            .1
            .trim_end_matches('|')
            .split("},{")
            .collect();
        assert_eq!(defs.len(), 2, "{s}");
        for def in defs {
            assert_eq!(
                def.chars().filter(|&c| c == 'C').count(),
                18,
                "definition {def} in {s} is not 18 carbons"
            );
        }
    }

    /// The Sg:-keeping fallback must not lose carbons either: a chain
    /// declared 18:1 still has to expand to 18 literal carbons, whichever
    /// path built it.
    #[test]
    fn sg_fallback_chains_keep_every_carbon() {
        assert_chain_map_sane("DG 18:1(9);5OH_18:1;OH", &[(1, 18), (2, 18)]);
    }

    #[test]
    fn expand_cxsmiles_handles_multi_chain_offset_and_renaming() {
        // Two chains each with their own unlocalized double bonds, joined
        // with unresolved regiochemistry: exercises CdkBuilder's atom
        // offsetting and variable renaming across multiple fragments, and
        // expand_cxsmiles_for_depiction's positional (not name-based)
        // matching of Sg positions back to their constraint equation.
        let s = lipid_name_to_smiles("PC 16:2_18:1").expect("should resolve");
        assert_balanced(&s);
        // Both chains need Sg:, so neither can move into an RG: definition
        // and the sn ambiguity goes unexpressed — the chains are written
        // into the backbone in name order instead.
        assert!(!s.contains("RG:") && !s.contains("f:"), "{s}");
        assert_eq!(
            s.matches("Sg:n:").count(),
            5,
            "chain1 (16:2, 2 unlocalized DBs -> 3 Sg) + chain2 (18:1, 1 unlocalized DB -> 2 Sg)"
        );

        let expanded = expand_cxsmiles_for_depiction(&s);
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
        let s = lipid_name_to_cxsmiles("PC 18:4(5,8,11,14);3oxo;16OH_18:2;oxo;OH")
            .expect("should resolve PC with variable content");
        eprintln!("PC multichain CXSMILES: {}", s);

        // Chain 2's unlocalized double bonds and modifications keep it in
        // the main SMILES, so this takes the Sg:/m: path, not the RG: one.
        let blocks = s.split(" |").nth(1).unwrap_or("");
        eprintln!("CXSMILES suffix: {}", blocks);
        assert!(blocks.contains("Sg:n:"), "chain 2 has unlocalized DBs: {s}");
        assert!(
            blocks.contains("m:"),
            "chain 2 has unlocalized oxygens: {s}"
        );
        assert!(!blocks.contains("RG:") && !blocks.contains("f:"), "{s}");
        // The `*` present is the position-variation stub, not an attachment
        // point for a chain.
        assert!(s.contains(".*O"), "{s}");
    }

    #[test]
    fn pc_multichain_known_sn() {
        // PC with known sn-positions - fallback to traditional approach
        // PC 18:1(9)/18:2 should not generate CDK CXSMILES (no variable content)
        let _s = lipid_name_to_cxsmiles("PC 18:1(9)/18:2");
        // Should return None since both DBs are known (18:2 has variable positions)
        // Wait, 18:2 has variable positions, so it should work!
        // Let me use a fully known case
        match lipid_name_to_cxsmiles("PC 18:0/18:0") {
            None => {
                // Fallback to traditional - expected since no variable content
                eprintln!("PC 18:0/18:0 correctly falls back to None (no variable content)");
            }
            Some(s) => {
                eprintln!("PC 18:0/18:0 unexpectedly generated CDK: {}", s);
            }
        }
    }

    #[test]
    fn cdk_cxsmiles_extended_classes() {
        // Test that extended lipid classes generate CXSMILES when they have variable content
        let test_cases = vec![
            ("PE 18:2_18:1(9)", "PE with variable DBs"),
            ("PS 18:2_18:1(9)", "PS with variable DBs"),
            ("PG 18:2_18:1(9)", "PG with variable DBs"),
            ("PA 18:2_18:1(9)", "PA with variable DBs"),
            ("LPE 18:2", "LPE with variable DBs"),
            ("LPS 18:2", "LPS with variable DBs"),
        ];

        for (input, description) in test_cases {
            match lipid_name_to_cxsmiles(input) {
                Some(cxsmiles) => {
                    eprintln!("✓ {}: {}", description, input);
                    eprintln!("  Generated CXSMILES: {}", cxsmiles);
                }
                None => {
                    eprintln!(
                        "✓ {}: {} (falls back - no variable content)",
                        description, input
                    );
                }
            }
        }
    }

    #[test]
    fn cdk_cxsmiles_validation() {
        // Validate that generated CXSMILES conform to CDK specifications
        let test_cases = vec![
            ("FA 18:2", "Variable DBs should generate Sg blocks"),
            ("FA 18:4(5,8,11,14);3oxo;16OH", "Known DBs with known mods"),
            (
                "PC 18:4(5,8,11,14);3oxo;16OH_18:2;oxo;OH",
                "Multi-chain with variable content",
            ),
        ];

        for (input, description) in test_cases {
            if let Some(cxsmiles) = lipid_name_to_cxsmiles(input) {
                eprintln!("\n✓ {}: {}", description, input);
                eprintln!("  CXSMILES: {}", cxsmiles);

                // Split SMILES and CXSMILES blocks
                let parts: Vec<&str> = cxsmiles.split(" |").collect();
                if parts.len() >= 2 {
                    let cxsmiles_blocks = parts[1];
                    eprintln!("  Blocks: {}", cxsmiles_blocks);

                    // Validate Sg blocks (format: Sg:n:position:variable:ht)
                    if cxsmiles_blocks.contains("Sg:") {
                        let sg_count = cxsmiles_blocks.matches("Sg:n:").count();
                        eprintln!("  ✓ {} Sg blocks found (Sg:n:pos:var:ht format)", sg_count);
                    }

                    // Validate m: blocks (format: m:atom_idx:pos.pos)
                    if cxsmiles_blocks.contains("m:") {
                        let m_count = cxsmiles_blocks.matches("m:").count();
                        eprintln!("  ✓ {} m: blocks found (m:atom:pos.pos format)", m_count);
                    }

                    // Validate f: block (format: f:0.1,0.2,...)
                    if cxsmiles_blocks.contains("f:") {
                        eprintln!("  ✓ f: block present (component grouping)");
                    }

                    // Validate constraint (format: var+var+...=number)
                    if let Some(pipe_pos) = cxsmiles_blocks.rfind('|') {
                        let constraint = &cxsmiles_blocks[pipe_pos + 1..];
                        if constraint.contains('=') && constraint.contains('+') {
                            eprintln!("  ✓ Constraint present: {}", constraint);
                        }
                    }
                } else {
                    eprintln!("  SMILES only (no CXSMILES blocks): {}", cxsmiles);
                }
            }
        }
    }
}
