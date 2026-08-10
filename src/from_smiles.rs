//! The reverse direction: a SMILES/CXSMILES string back to a Shorthand2020
//! name.
//!
//! # What this inverts, and what it does not
//!
//! This reads back what [`crate::name2smiles`] writes. It is not a general
//! SMILES parser and does not try to be one: it recognizes this generator's
//! own headgroup templates and chain layout, and returns `None` for anything
//! else. Ingesting arbitrary third-party structures needs real substructure
//! perception, which a zero-dependency string manipulator has no business
//! pretending to do.
//!
//! # Every answer is proved before it is returned
//!
//! Reading a structure back is exactly the kind of index arithmetic that goes
//! subtly wrong and stays wrong, so nothing here is trusted on its own. The
//! name this module assembles is fed back through [`lipid_name_to_smiles`]
//! and compared to the input; unless it regenerates the string *exactly*, the
//! answer is discarded and `None` is returned.
//!
//! That makes the function self-checking in the strongest available sense: a
//! bug in this file can cost coverage, but it cannot produce a name that
//! means something other than the structure handed in. It also lets the
//! layout table below be tried speculatively — a template that does not
//! belong simply fails to verify.
//!
//! # Where the round trip is lossy
//!
//! [`crate::name2smiles`] is not injective, so `name -> SMILES -> name` does
//! not always land on the name it started from. It always lands on a name
//! that regenerates the same string, which is the property worth having:
//!
//! * A trailing empty chain slot is dropped, so `DG 16:0/0:0` and `MG 16:0`
//!   produce the same string and both come back as `MG 16:0`.
//! * Older spellings canonicalize: `;ep(5)` comes back as `;5Ep`, `;OH(3)` as
//!   `;3OH`.

use crate::smiles::{
    count_atoms, lipid_name_to_smiles, trailer_equations, trailer_is_swappable, SUBSTITUENTS,
};

/// A structure read back far enough to name one chain.
#[derive(Debug, Default, Clone)]
struct Chain {
    /// `""`, `"O-"`, `"P-"`, `"d"` or `"t"`.
    prefix: &'static str,
    carbon: u32,
    /// Localized double bonds: Δ-position and geometry, if the string
    /// carried one.
    db: Vec<(u32, Option<char>)>,
    /// Double bonds the `Sg:` scaffold stands in for — counted, never placed.
    unlocalized_db: usize,
    /// `(position, abbreviation)`, position `0` meaning an `m:` block said
    /// "present, not determined".
    mods: Vec<(u32, &'static str)>,
    /// `(first carbon, last carbon, is_epoxide)`.
    rings: Vec<(u32, u32, bool)>,
}

/// The CXSMILES tail, split into the parts that mean something here: the
/// standard `|...|` blocks, and the trailing token list that is ours
/// (`dev/extension.md`).
struct Cx<'a> {
    base: &'a str,
    /// Atom index of each `Sg:n:` marker, in emission order.
    sg: Vec<usize>,
    /// `(floating atom index, candidate atom indices)` for each `m:` block.
    m: Vec<(usize, Vec<usize>)>,
    /// Whether a `swappable(...)` token says the sn assignment shown is one
    /// arbitrary choice — i.e. whether the name joined its chains with `_`.
    swappable: bool,
    /// `(term count, sum)` per `constrain(...)` token, in order. They pair
    /// with `sg` positionally — see `expand_cxsmiles_for_depiction`.
    equations: Vec<(usize, usize)>,
}

impl<'a> Cx<'a> {
    fn split(smi: &'a str) -> Cx<'a> {
        let Some((base, tail)) = smi.split_once(" |") else {
            return Cx {
                base: smi,
                sg: Vec::new(),
                m: Vec::new(),
                swappable: false,
                equations: Vec::new(),
            };
        };
        let (blocks, trailer) = match tail.split_once('|') {
            Some((b, c)) => (b, c),
            None => (tail, ""),
        };

        let mut sg = Vec::new();
        let mut m = Vec::new();
        for block in blocks.split(',') {
            if let Some(rest) = block.strip_prefix("Sg:n:") {
                if let Some(Ok(atom)) = rest.split(':').next().map(str::parse) {
                    sg.push(atom);
                }
            } else if let Some(rest) = block.strip_prefix("m:") {
                if let Some((idx, sites)) = rest.split_once(':') {
                    if let Ok(idx) = idx.parse() {
                        m.push((
                            idx,
                            sites.split('.').filter_map(|s| s.parse().ok()).collect(),
                        ));
                    }
                }
            }
        }

        // The trailer grammar lives with the emitter, in `smiles.rs`.
        let swappable = trailer_is_swappable(trailer);
        let equations = trailer_equations(trailer);

        Cx {
            base,
            sg,
            m,
            swappable,
            equations,
        }
    }

    /// The `Sg:` scaffold belonging to the chain occupying atoms
    /// `span`: how many double bonds it stands in for, and how many carbons
    /// expansion would add. Equations pair with markers in emission order,
    /// so this walks both together.
    fn scaffold(&self, span: std::ops::Range<usize>) -> (usize, u32) {
        let mut consumed = 0;
        for &(terms, sum) in &self.equations {
            if consumed + terms > self.sg.len() {
                break;
            }
            let markers = &self.sg[consumed..consumed + terms];
            consumed += terms;
            if markers.iter().all(|a| span.contains(a)) {
                // The scaffold writes one marker per double bond plus one,
                // and every variable is at least 1, so the carbons it hides
                // are the constraint sum less the marker count.
                return (terms - 1, (sum - terms) as u32);
            }
        }
        (0, 0)
    }
}

/// Reads one chain fragment — the SMILES written for a single acyl/alkyl
/// chain, plus any `.*X` components its `m:` blocks point at — into the chain
/// it was built from. `atom_base` is the fragment's first atom index in the
/// whole string, which is what the `Sg:`/`m:` blocks are numbered against.
///
/// `first_carbon` is the chain position the fragment's first carbon holds: 1
/// for an ordinary chain, 3 for a sphingoid base's tail, whose C1 and C2 live
/// in the headgroup template.
fn read_chain(text: &str, atom_base: usize, first_carbon: u32, cx: &Cx) -> Option<Chain> {
    let mut parts = text.split('.');
    let main = parts.next()?;

    let mut chain = Chain::default();
    let chars: Vec<char> = main.chars().collect();
    let mut i = 0;
    let mut atom = atom_base;
    let mut carbon = first_carbon;
    // Bond text sitting between the previous carbon and this one, and the
    // position of that previous carbon.
    let mut pending_bond = String::new();
    let mut labels: Vec<(usize, u32)> = Vec::new();
    let mut bonds: Vec<(u32, String)> = Vec::new();
    let mut carbon_atoms: Vec<usize> = Vec::new();

    while i < chars.len() {
        match chars[i] {
            'C' => {
                if carbon > first_carbon {
                    bonds.push((carbon - 1, std::mem::take(&mut pending_bond)));
                }
                carbon_atoms.push(atom);
                atom += 1;
                i += 1;
                // Ring-closure labels bind to the atom they follow.
                while i < chars.len() && chars[i] == '%' {
                    let start = i + 1;
                    let mut end = start;
                    while end < chars.len() && chars[end].is_ascii_digit() {
                        end += 1;
                    }
                    labels.push((
                        chars[start..end].iter().collect::<String>().parse().ok()?,
                        carbon,
                    ));
                    i = end;
                }
                // Then its substituent branches.
                while i < chars.len() && chars[i] == '(' {
                    let (branch, next) = read_branch(&chars, i)?;
                    atom += count_atoms(&branch);
                    if carbon == 1 && branch == "=O" && first_carbon == 1 {
                        // C1's carbonyl is the ester/amide linkage itself,
                        // not an `oxo` group written on the chain.
                        chain.prefix = "";
                    } else {
                        chain.mods.push((carbon, abbreviation(&branch)?));
                    }
                    i = next;
                }
                carbon += 1;
            }
            'O' => {
                // The only bare heteroatom inside a chain fragment is an
                // epoxide's bridging oxygen.
                pending_bond.push('O');
                atom += 1;
                i += 1;
            }
            c @ ('=' | '/' | '\\' | '#') => {
                pending_bond.push(c);
                i += 1;
            }
            _ => return None,
        }
    }
    chain.carbon = carbon - 1;
    if chain.carbon == 0 {
        return None;
    }

    // An acyl chain is the one whose C1 carries the carbonyl; anything else
    // is an ether linkage, vinyl if its first bond is the mandatory C1=C2.
    let acyl = main.starts_with("C(=O)") && first_carbon == 1;
    if first_carbon == 1 && !acyl {
        chain.prefix = if bonds.first().is_some_and(|(_, b)| b.contains('=')) {
            "P-"
        } else {
            "O-"
        };
    }

    for (k, bond) in &bonds {
        if bond.contains('=') {
            chain.db.push((*k, None));
        }
    }
    // Geometry is written as a marker on the bond *after* the double bond:
    // `\` closes a cis pair, `/` a trans one.
    for (pos, geom) in chain.db.iter_mut() {
        if let Some((_, after)) = bonds.iter().find(|(k, _)| *k == *pos + 1) {
            *geom = match () {
                _ if after.contains('\\') => Some('Z'),
                _ if after.contains('/') => Some('E'),
                _ => None,
            };
        }
    }

    // Ring-closure labels used exactly twice are a ring; an epoxide is the
    // adjacent-carbon case whose bond text carries the bridging oxygen.
    labels.sort();
    for pair in labels.chunks(2) {
        let [(la, a), (lb, b)] = pair else {
            return None;
        };
        if la != lb {
            return None;
        }
        let bridged = bonds
            .iter()
            .any(|(k, bond)| k == a.min(b) && bond.contains('O'));
        chain.rings.push((*a.min(b), *a.max(b), bridged));
    }
    // A ring bond is not a chain double bond the name should place.
    chain.db.retain(|(pos, _)| {
        !chain
            .rings
            .iter()
            .any(|(s, e, bridged)| *bridged && pos >= s && pos < e)
    });

    let span = atom_base..atom;
    let (unlocalized, hidden) = cx.scaffold(span.clone());
    chain.unlocalized_db = unlocalized;
    chain.carbon += hidden;
    // The scaffold's own double bonds are the unlocalized ones, and they are
    // always the last written, so they are dropped from the placed list.
    for _ in 0..unlocalized {
        chain.db.pop();
    }

    // Each `.*X` component is a group the name declared without a position.
    // Its `m:` block has to point at this chain's carbons for it to be ours.
    let mut floating = atom;
    for part in parts {
        let branch = part.strip_prefix('*')?;
        let mine = cx
            .m
            .iter()
            .any(|(idx, sites)| *idx == floating && sites.iter().all(|s| carbon_atoms.contains(s)));
        if !mine {
            return None;
        }
        chain.mods.push((0, abbreviation(branch)?));
        floating += count_atoms(part);
    }

    chain.db.sort();
    chain.mods.sort();
    Some(chain)
}

/// Reads the parenthesized branch starting at `chars[open]`, returning its
/// contents and the index just past its closing paren.
fn read_branch(chars: &[char], open: usize) -> Option<(String, usize)> {
    let mut depth = 0;
    for (i, &c) in chars.iter().enumerate().skip(open) {
        match c {
            '(' => depth += 1,
            ')' => {
                depth -= 1;
                if depth == 0 {
                    return Some((chars[open + 1..i].iter().collect(), i + 1));
                }
            }
            _ => {}
        }
    }
    None
}

/// The Table 1A abbreviation whose branch is `branch`.
fn abbreviation(branch: &str) -> Option<&'static str> {
    SUBSTITUENTS
        .iter()
        .find(|(_, b)| *b == branch)
        .map(|(abbr, _)| *abbr)
}

/// Writes a chain back out as its shorthand token, e.g. `18:1(9Z);5OH`.
fn format_chain(chain: &Chain) -> String {
    // A ring's double bonds are declared inside the ring and counted there,
    // so they belong to neither the chain's `C:DBE` count nor its position
    // list. A plasmenyl chain's C1=C2 is likewise implied by the `P-` prefix.
    let in_ring = |pos: u32| {
        chain
            .rings
            .iter()
            .any(|(s, e, bridged)| !*bridged && pos >= *s && pos < *e)
    };
    let placed: Vec<&(u32, Option<char>)> = chain
        .db
        .iter()
        .filter(|(p, _)| !(in_ring(*p) || (chain.prefix == "P-" && *p == 1)))
        .collect();

    let mut out = format!(
        "{}{}:{}",
        chain.prefix,
        chain.carbon,
        placed.len() + chain.unlocalized_db
    );
    if !placed.is_empty() {
        let list: Vec<String> = placed
            .iter()
            .map(|(p, g)| match g {
                Some(g) => format!("{p}{g}"),
                None => p.to_string(),
            })
            .collect();
        out.push_str(&format!("({})", list.join(",")));
    }

    // Rings come before the other groups, with any group inside their span
    // written within their brackets.
    for &(start, end, bridged) in &chain.rings {
        if bridged {
            out.push_str(&format!(";{start}Ep"));
            continue;
        }
        let ring_db: Vec<String> = chain
            .db
            .iter()
            .filter(|(p, _)| *p >= start && *p < end)
            .map(|(p, g)| match g {
                Some(g) => format!("{p}{g}"),
                None => p.to_string(),
            })
            .collect();
        let mut inner = format!("{start}-{end}cy{}:{}", end - start + 1, ring_db.len());
        if !ring_db.is_empty() {
            inner.push_str(&format!("({})", ring_db.join(",")));
        }
        for (pos, abbr) in chain.mods.iter().filter(|(p, _)| *p >= start && *p <= end) {
            inner.push_str(&format!(";{pos}{abbr}"));
        }
        out.push_str(&format!(";[{inner}]"));
    }

    // Then the rest, grouped by abbreviation in the table's own order so a
    // name comes back out the way the paper writes it.
    for (abbr, _) in SUBSTITUENTS {
        let positioned: Vec<String> = chain
            .mods
            .iter()
            .filter(|(p, a)| a == abbr && *p != 0 && !in_ring_span(chain, *p))
            .map(|(p, a)| format!("{p}{a}"))
            .collect();
        if !positioned.is_empty() {
            out.push_str(&format!(";{}", positioned.join(",")));
        }
        let loose = chain
            .mods
            .iter()
            .filter(|(p, a)| a == abbr && *p == 0)
            .count();
        if loose > 0 {
            // The paper parenthesizes a group when it occurs more than once,
            // and whenever its abbreviation carries a digit of its own.
            if loose > 1 || abbr.contains(|c: char| c.is_ascii_digit()) {
                out.push_str(&format!(";({abbr})"));
                if loose > 1 {
                    out.push_str(&loose.to_string());
                }
            } else {
                out.push_str(&format!(";{abbr}"));
            }
        }
    }
    out
}

fn in_ring_span(chain: &Chain, pos: u32) -> bool {
    chain
        .rings
        .iter()
        .any(|(s, e, bridged)| !*bridged && pos >= *s && pos <= *e)
}

/// One class's resolved template: the literal text the builder writes, split
/// at its chain slots, and the sn position each gap holds.
///
/// A slot's leading ester/ether `O` belongs to the segment before it, so an
/// unfilled position is simply an empty gap — which is exactly the bare `O`
/// `slot_fragment_cdk` writes for it.
struct Layout {
    class: &'static str,
    segments: &'static [&'static str],
    sn: &'static [usize],
}

/// Every glycerol-backboned class shares one shape; which of MG/DG/TG it is
/// follows from how far along the last filled slot sits.
const GLYCEROL: Layout = Layout {
    class: "",
    segments: &["C(CO", ")(O", ")CO", ""],
    sn: &[3, 2, 1],
};

const GPL_TAILS: &[(&str, &str)] = &[
    ("PC", "OP(=O)([O-])OCC[N+](C)(C)C"),
    ("PE", "OP(=O)(O)OCCN"),
    ("PS", "OP(=O)(O)OCC(N)C(=O)O"),
    ("PG", "OP(=O)(O)OCC(O)CO"),
    ("PI", "OP(=O)(O)OC1C(O)C(O)C(O)C(O)C1O"),
    ("PA", "OP(=O)(O)O"),
];

/// `(class, headgroup tail, carries an N-acyl chain)`.
const SPHINGO: &[(&str, &str, bool)] = &[
    ("SM", "OP(=O)([O-])OCC[N+](C)(C)C", true),
    ("HexCer", "OC1OC(CO)C(O)C(O)C1O", true),
    ("IPC", "OP(=O)(O)OC1C(O)C(O)C(O)C(O)C1O", true),
    ("CerP", "OP(=O)(O)O", true),
    ("Cer", "O", true),
    ("S1P", "OP(=O)(O)O", false),
    ("Sph", "O", false),
];

const SINGLE_CHAIN: &[Layout] = &[
    Layout {
        class: "AMP-FA",
        segments: &["[n+]1ccccc1-c1ccc(CN", ")cc1"],
        sn: &[1],
    },
    Layout {
        class: "CE",
        segments: &["C12(CC=C3CC(O", ")CCC3(C)C1CCC1(C)C(C(C)CCCC(C)C)CCC21)"],
        sn: &[1],
    },
    Layout {
        class: "CAR",
        segments: &["O(", ")[C@H](CC(=O)[O-])C[N+](C)(C)C"],
        sn: &[1],
    },
    Layout {
        class: "NAE",
        segments: &["OCCN", ""],
        sn: &[1],
    },
    Layout {
        class: "FA",
        segments: &["O", ""],
        sn: &[1],
    },
];

/// Splits `base` along a layout's literal segments, returning the text and
/// starting atom index of each chain gap.
///
/// A gap ends at the first point where the remaining text begins with the
/// next segment *and* the parentheses opened inside the gap have all closed —
/// a chain's own `(=O)` branches close before the backbone's do.
fn split_layout<'a>(base: &'a str, segments: &[&str]) -> Option<Vec<(&'a str, usize)>> {
    let mut gaps = Vec::new();
    let mut rest = base;
    let mut atom = 0;

    let first = segments.first()?;
    rest = rest.strip_prefix(first)?;
    atom += count_atoms(first);

    for next in &segments[1..] {
        let end = if next.is_empty() {
            // The last segment: the gap runs to the end of the string.
            rest.len()
        } else if rest.starts_with(*next) {
            // The next segment starts immediately — an unfilled position.
            0
        } else {
            let mut depth = 0i32;
            let mut found = None;
            for (i, c) in rest.char_indices() {
                match c {
                    '(' => depth += 1,
                    ')' => depth -= 1,
                    _ => {}
                }
                if depth < 0 {
                    break;
                }
                let at = i + c.len_utf8();
                if depth == 0 && rest[at..].starts_with(*next) {
                    found = Some(at);
                    break;
                }
            }
            found?
        };
        gaps.push((&rest[..end], atom));
        atom += count_atoms(&rest[..end]);
        rest = &rest[end..];
        rest = rest.strip_prefix(*next)?;
        atom += count_atoms(next);
    }
    if rest.is_empty() {
        Some(gaps)
    } else {
        None
    }
}

/// Assembles the chain tokens into sn order, writing `0:0` for a gap that
/// some later position fills. The separator is `_` when a `swappable` token
/// said the sn assignment shown is arbitrary, `/` when it is determined.
fn join_slots(mut slots: Vec<(usize, Option<String>)>, sep: &str) -> Option<(usize, String)> {
    slots.sort_by_key(|(sn, _)| *sn);
    let last = slots.iter().rposition(|(_, c)| c.is_some())?;
    let tokens: Vec<String> = slots[..=last]
        .iter()
        .map(|(_, c)| c.clone().unwrap_or_else(|| "0:0".to_string()))
        .collect();
    Some((last + 1, tokens.join(sep)))
}

/// Reads a SMILES/CXSMILES string back into a Shorthand2020 name, or `None`
/// when no name this generator would accept regenerates it exactly.
pub fn smiles_to_name(smi: &str) -> Option<String> {
    let smi = smi.trim();
    let cx = Cx::split(smi);

    candidates(&cx)
        .into_iter()
        .find(|candidate| lipid_name_to_smiles(candidate).as_deref() == Some(smi))
}

/// Every name worth testing against the input, cheapest and most specific
/// first. Each is only a hypothesis — `smiles_to_name` keeps the one that
/// regenerates the string.
fn candidates(cx: &Cx) -> Vec<String> {
    let mut out = Vec::new();
    let sep = if cx.swappable { "_" } else { "/" };

    // Sterol has no chains at all.
    out.push("ST".to_string());

    let read = |text: &str, atom: usize| -> Option<String> {
        (!text.is_empty()).then(|| read_chain(text, atom, 1, cx).as_ref().map(format_chain))?
    };

    for layout in SINGLE_CHAIN {
        if let Some(gaps) = split_layout(cx.base, layout.segments) {
            if let Some(token) = read(gaps[0].0, gaps[0].1) {
                out.push(format!("{} {token}", layout.class));
            }
        }
    }

    if let Some(gaps) = split_layout(cx.base, GLYCEROL.segments) {
        let slots: Vec<(usize, Option<String>)> = GLYCEROL
            .sn
            .iter()
            .zip(&gaps)
            .map(|(sn, (text, atom))| (*sn, read(text, *atom)))
            .collect();
        if let Some((filled, tokens)) = join_slots(slots, sep) {
            out.push(format!("{} {tokens}", ["MG", "DG", "TG"][filled - 1]));
        }
    }

    for (class, tail) in GPL_TAILS {
        let segments = [format!("C(C{tail})(O"), ")CO".to_string(), String::new()];
        let segments: Vec<&str> = segments.iter().map(String::as_str).collect();
        if let Some(gaps) = split_layout(cx.base, &segments) {
            let slots: Vec<(usize, Option<String>)> = [2, 1]
                .iter()
                .zip(&gaps)
                .map(|(sn, (text, atom))| (*sn, read(text, *atom)))
                .collect();
            if let Some((filled, tokens)) = join_slots(slots, sep) {
                // One chain at sn1 and nothing at sn2 is the lyso form.
                out.push(format!(
                    "{}{class} {tokens}",
                    if filled == 1 { "L" } else { "" }
                ));
            }
        }
    }

    let cl = ["C(COP(=O)(O)OCC(O", ")CO", ")(O)COP(=O)(O)OCC(O", ")CO", ""];
    if let Some(gaps) = split_layout(cx.base, &cl) {
        let slots: Vec<(usize, Option<String>)> = [2, 1, 4, 3]
            .iter()
            .zip(&gaps)
            .map(|(sn, (text, atom))| (*sn, read(text, *atom)))
            .collect();
        if let Some((_, tokens)) = join_slots(slots, sep) {
            out.push(format!("CL {tokens}"));
        }
    }

    out.extend(sphingoid_candidates(cx));
    out
}

/// Sphingoids need their own reader: the base's C1 and C2 are written by the
/// template, its C3 (and C4, on a triol) carry hydroxyls the `d`/`t` prefix
/// implies rather than states, and its tail picks up at C3.
fn sphingoid_candidates(cx: &Cx) -> Vec<String> {
    let mut out = Vec::new();
    for (class, tail, has_n_acyl) in SPHINGO {
        let head = format!("C(C{tail})(N");
        let Some(rest) = cx.base.strip_prefix(&head) else {
            continue;
        };
        let mut atom = count_atoms(&head);

        // The N-acyl chain, if this class has one, runs to the paren that
        // closes the branch the template opened.
        let mut depth = 1i32;
        let mut split = None;
        for (i, c) in rest.char_indices() {
            match c {
                '(' => depth += 1,
                ')' => {
                    depth -= 1;
                    if depth == 0 {
                        split = Some(i);
                        break;
                    }
                }
                _ => {}
            }
        }
        let Some(split) = split else { continue };
        let (acyl_text, base_text) = (&rest[..split], &rest[split + 1..]);
        if acyl_text.is_empty() == *has_n_acyl {
            continue;
        }

        let acyl = if *has_n_acyl {
            let Some(c) = read_chain(acyl_text, atom, 1, cx) else {
                continue;
            };
            atom += count_atoms(acyl_text);
            Some(format_chain(&c))
        } else {
            None
        };

        let Some(mut base) = read_chain(base_text, atom, 3, cx) else {
            continue;
        };
        // The tail was read with C3 as its first carbon, so `carbon` is
        // already the base's full length — C1 and C2 live in the template.
        // The hydroxyls at C3 (and C4 for a triol) are what the prefix
        // letter says rather than something the name spells out.
        let triol = base.mods.contains(&(4, "OH"));
        if !base.mods.contains(&(3, "OH")) {
            continue;
        }
        base.prefix = if triol { "t" } else { "d" };
        base.mods.retain(|m| *m != (3, "OH") && *m != (4, "OH"));

        let token = format_chain(&base);
        out.push(match acyl {
            Some(acyl) => format!("{class} {token}/{acyl}"),
            None => format!("{class} {token}"),
        });
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn round(name: &str) -> Option<String> {
        smiles_to_name(&lipid_name_to_smiles(name)?)
    }

    #[test]
    fn plain_chains_come_back_exactly() {
        for name in [
            "FA 18:0",
            "FA 18:1(9)",
            "FA 18:1(9Z)",
            "FA 18:2(9Z,12Z)",
            "PC 16:0/18:1(9Z)",
            "TG 16:0/18:1(9)/18:2(9,12)",
            "CL 16:0/18:1(9)/16:0/18:1(9)",
            "CE 18:1(9Z)",
            "CAR 18:1(9)",
            "NAE 20:4(5,8,11,14)",
            "AMP-FA 20:4(5,8,11,14);15OH",
            "MG 18:1(9)",
            "LPC 18:1(9)",
            "PC O-16:0/18:1(9)",
            "PC P-16:0/18:1(9)",
            "ST",
        ] {
            assert_eq!(round(name).as_deref(), Some(name), "{name}");
        }
    }

    /// The whole point of the exercise: what was *not* determined has to
    /// survive the trip as undetermined, rather than coming back as a
    /// position someone will later read as a measurement.
    #[test]
    fn undetermined_things_stay_undetermined() {
        for name in [
            "FA 18:1",               // unlocalized double bond -> Sg:
            "FA 20:4",               // four of them
            "FA 18:2(9)",            // one placed, one not
            "FA 18:0;OH",            // unlocalized modification -> m:
            "FA 18:1;OH",            // both at once
            "FA 20:3(5,8,11);(OH)2", // two unlocalized groups
            "PC 16:0_18:1(9)",       // unresolved sn -> RG:
            "TG 16:0_18:1(9)_18:0",  // three-way
            "CL 16:0_18:1(9)_16:0_18:1(9)",
        ] {
            assert_eq!(round(name).as_deref(), Some(name), "{name}");
        }
    }

    #[test]
    fn sphingoid_bases_recover_their_prefix_and_n_acyl() {
        for name in [
            "Cer d18:1(4)/16:0",
            "SM d18:1(4)/16:0",
            "HexCer d18:1(4)/16:0",
            "Sph d18:1(4)",
            "S1P d18:1(4)",
            "Cer d18:1(4)/16:0;9Ep",
        ] {
            assert_eq!(round(name).as_deref(), Some(name), "{name}");
        }
    }

    #[test]
    fn table_1a_groups_and_table_1b_rings_survive() {
        for name in [
            "FA 18:0;5Me",
            "FA 16:0;3Me,7Me,11Me,15Me",
            "FA 18:0;5Br",
            "FA 18:0;5CN",
            "FA 18:1(9);12NO2",
            "FA 18:1(9);(NO2)",
            "FA 20:3(5Z,13E);11OH,15OH;9oxo",
            "FA 18:0;9Ep",
            "FA 19:0;[11-13cy3:0]",
            "FA 19:0;[9-11cy3:1(9)]",
            "FA 18:1(6Z);[14-18cy5:1(15)]",
        ] {
            assert_eq!(round(name).as_deref(), Some(name), "{name}");
        }
    }

    /// Not a general SMILES parser, and honest about it.
    #[test]
    fn foreign_structures_are_refused() {
        for smi in [
            "",
            "c1ccccc1",                                 // benzene
            "CC(=O)O",                                  // acetic acid, no headgroup template
            "OC(=O)CCCCCCCCCCCCCCCCC.O",                // a stearate with something extra
            "C(COP(=O)([O-])OCC[N+](C)(C)C)(OC(=O)CCC", // truncated
            "OC(=O)CCCC[Se]CCCC",                       // element this crate never writes
        ] {
            assert_eq!(smiles_to_name(smi), None, "{smi:?} should be refused");
        }
    }

    /// A returned name must regenerate the exact string it was read from.
    /// Tampering with a structure therefore has only two honest outcomes: a
    /// name for the *tampered* structure, or nothing at all — never the name
    /// of the structure it was derived from.
    #[test]
    fn a_tampered_structure_never_returns_the_original_name() {
        let name = "PC 16:0/18:1(9)";
        let original = lipid_name_to_smiles(name).unwrap();

        for tampered in [
            original.replacen("CCCCCCCC=C", "CCCCCCC=C", 1), // move the double bond
            original.replacen("COC(=O)CCCC", "COC(=O)CCC", 1), // shorten a chain
            original.replace("[N+](C)(C)C", "[N+](C)(C)CC"), // alter the headgroup
        ] {
            assert_ne!(tampered, original, "tamper should have changed something");
            match smiles_to_name(&tampered) {
                None => {}
                Some(recovered) => {
                    assert_ne!(recovered, name, "{tampered}: returned the untampered name");
                    assert_eq!(
                        lipid_name_to_smiles(&recovered).as_deref(),
                        Some(tampered.as_str()),
                        "{tampered}: recovered name does not regenerate it"
                    );
                }
            }
        }
    }

    /// A chain read back and written out again must be a fixed point: if the
    /// first trip canonicalized anything, the second must change nothing.
    #[test]
    fn the_round_trip_is_idempotent() {
        for name in [
            "FA 18:0;ep(5)",
            "FA 18:0;OH(5)",
            "PC 16:0_18:1",
            "DG 16:0/0:0",
        ] {
            let once = round(name).unwrap_or_else(|| panic!("{name} should resolve"));
            let twice = round(&once).unwrap_or_else(|| panic!("{once} should resolve"));
            assert_eq!(once, twice, "{name} did not settle after one trip");
        }
    }
}
