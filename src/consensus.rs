//! The bracketed consensus tail, carried as trailer tokens.
//!
//! The EAD engines emit names with a per-position consensus in square
//! brackets:
//!
//! ```text
//! FA 20:4;11OH [DB sn1: Δ5 100%, Δ8 100%, Δ12 100% | Δ14 50% | Δ15 50%]
//! ```
//!
//! No structure format can hold a weighted call, so this used to be stripped
//! and thrown away. It is *metadata*, though, and the trailer after the
//! closing pipe is exactly where this crate's metadata lives — so it can be
//! carried there without any of it touching the CXSMILES.
//!
//! # What the tokens say
//!
//! | token | says |
//! |---|---|
//! | `dbPos(sn1:5@100,8@100)` | double bonds at Δ5 and Δ8, each called with 100% confidence |
//! | `dbPos(sn1:14@50\|15@50)` | one double bond, Δ14 or Δ15, evenly split |
//! | `mPos(sn1:11OH@100)` | the hydroxyl at position 11, called with 100% confidence |
//! | `mPos(OH1:11OH@50,13OH@50)` | the group on the stub labelled `OH1` is at 11 or 13 |
//!
//! `|` separates mutually exclusive alternatives, which is what it means in
//! the source tail: one feature, several candidate positions.
//!
//! # Two rules this obeys
//!
//! **It never contradicts the structure, only refines it.** A `dbPos` entry
//! for a position the SMILES already commits to records how sure that call
//! was. A set of `|` alternatives corresponds to a double bond the structure
//! left inside an `Sg:` run, or a group it left on an `m:` stub — the token
//! narrows "somewhere in this stretch" to "one of these, with these odds",
//! which is strictly more information and still not a determination.
//!
//! **It names things, never positions.** The anchor before the `:` is an
//! atom label the `|...|` block carries — `snN` for a chain, or the stub's
//! own label. An atom *index* would rot the first time a toolkit renumbered
//! the molecule, because nothing rewrites the title field. Verified: a label
//! on a floating `*` stub survives canonical rewriting still attached to its
//! own atom.
//!
//! # Syntax
//!
//! The source tail is written with spaces and `Δ` signs; the tokens are not.
//! The trailer must stay free of whitespace because `.smi` readers split the
//! line on it and keep field 1 as the ID, so `Δ5 100%` becomes `5@100`. The
//! tail is reconstructed in its original spelling on the way back out.

use std::fmt::Write as _;

/// One called position and the confidence behind it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Call {
    pub(crate) pos: u32,
    /// Percent, `0..=100`. The source writes `92%`; this holds `92`.
    pub(crate) percent: u32,
}

/// One feature's consensus on one chain.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Consensus {
    /// `"DB"` for a double bond, otherwise a Table 1A abbreviation (`"OH"`).
    pub(crate) kind: String,
    /// The sn position the feature sits on, from `sn1`, `sn2`, …
    pub(crate) sn: usize,
    /// `|`-separated groups of mutually exclusive candidates. A group with
    /// one entry is a call the engine was sure of; several entries are
    /// alternatives for a single feature.
    pub(crate) alternatives: Vec<Vec<Call>>,
}

impl Consensus {
    /// Whether this records a single determined position rather than a set of
    /// competing ones.
    pub(crate) fn is_localized(&self) -> bool {
        self.alternatives.iter().all(|group| group.len() == 1)
    }
}

/// Parses the bracketed tail into one entry per `;`-separated section.
///
/// Returns `None` for anything that is not a well-formed tail, so a name that
/// merely ends in a bracket is left alone.
pub(crate) fn parse_tail(tail: &str) -> Option<Vec<Consensus>> {
    let inner = tail.trim().strip_prefix('[')?.strip_suffix(']')?;
    let mut out = Vec::new();
    for section in inner.split(';') {
        let section = section.trim();
        if section.is_empty() {
            continue;
        }
        let (head, body) = section.split_once(':')?;
        let mut head = head.split_whitespace();
        let kind = head.next()?.to_string();
        let sn = head.next()?.strip_prefix("sn")?.parse().ok()?;
        if head.next().is_some() {
            return None;
        }

        let mut alternatives = Vec::new();
        for group in body.split('|') {
            let calls: Option<Vec<Call>> = group.split(',').map(parse_call).collect();
            let calls = calls?;
            if calls.is_empty() {
                return None;
            }
            alternatives.push(calls);
        }
        if alternatives.is_empty() {
            return None;
        }
        out.push(Consensus {
            kind,
            sn,
            alternatives,
        });
    }
    (!out.is_empty()).then_some(out)
}

/// `Δ5 100%` or `11 100%` into a [`Call`].
fn parse_call(text: &str) -> Option<Call> {
    let text = text.trim().trim_start_matches('\u{0394}');
    let (pos, percent) = text.split_once(char::is_whitespace)?;
    Some(Call {
        pos: pos.trim().parse().ok()?,
        percent: percent.trim().trim_end_matches('%').trim().parse().ok()?,
    })
}

/// The trailer tokens for a set of consensus entries.
///
/// `stub_label` is asked for the label of the `m:` stub belonging to a given
/// unlocalized group, so a token about a floating substituent anchors to that
/// substituent rather than merely to its chain.
pub(crate) fn tokens(
    entries: &[Consensus],
    mut stub_label: impl FnMut(&Consensus) -> Option<String>,
) -> Vec<String> {
    entries
        .iter()
        .map(|entry| {
            let anchor = stub_label(entry).unwrap_or_else(|| format!("sn{}", entry.sn));
            let body = entry
                .alternatives
                .iter()
                .map(|group| {
                    group
                        .iter()
                        .map(|call| {
                            let mut out = call.pos.to_string();
                            if entry.kind != "DB" {
                                out.push_str(&entry.kind);
                            }
                            let _ = write!(out, "@{}", call.percent);
                            out
                        })
                        .collect::<Vec<_>>()
                        .join(",")
                })
                .collect::<Vec<_>>()
                .join("|");
            let name = if entry.kind == "DB" { "dbPos" } else { "mPos" };
            format!("{name}({anchor}:{body})")
        })
        .collect()
}

/// Reads `dbPos`/`mPos` tokens back into consensus entries.
///
/// `sn_of` resolves an anchor label to the sn position it belongs to, so a
/// stub-anchored token recovers the chain it was written for.
pub(crate) fn from_tokens(
    trailer: &str,
    mut sn_of: impl FnMut(&str) -> Option<usize>,
) -> Vec<Consensus> {
    let mut out = Vec::new();
    for token in trailer.split(';') {
        let token = token.trim();
        let Some((name, rest)) = token.split_once('(') else {
            continue;
        };
        if name != "dbPos" && name != "mPos" {
            continue;
        }
        let Some((anchor, body)) = rest.trim_end_matches(')').split_once(':') else {
            continue;
        };
        let Some(sn) = sn_of(anchor) else { continue };

        let mut kind = String::from("DB");
        let mut alternatives = Vec::new();
        let mut ok = true;
        for group in body.split('|') {
            let mut calls = Vec::new();
            for item in group.split(',') {
                let Some((head, percent)) = item.split_once('@') else {
                    ok = false;
                    break;
                };
                let digits = head.trim_end_matches(|c: char| c.is_ascii_alphabetic());
                if name == "mPos" {
                    kind = head[digits.len()..].to_string();
                }
                match (digits.parse(), percent.parse()) {
                    (Ok(pos), Ok(percent)) => calls.push(Call { pos, percent }),
                    _ => {
                        ok = false;
                        break;
                    }
                }
            }
            if !ok || calls.is_empty() {
                ok = false;
                break;
            }
            alternatives.push(calls);
        }
        if ok {
            out.push(Consensus {
                kind,
                sn,
                alternatives,
            });
        }
    }
    out
}

/// Writes consensus entries back as the bracketed tail they came from.
pub(crate) fn format_tail(entries: &[Consensus]) -> Option<String> {
    if entries.is_empty() {
        return None;
    }
    let sections: Vec<String> = entries
        .iter()
        .map(|entry| {
            let groups: Vec<String> = entry
                .alternatives
                .iter()
                .map(|group| {
                    group
                        .iter()
                        .map(|call| {
                            let delta = if entry.kind == "DB" { "\u{0394}" } else { "" };
                            format!("{delta}{} {}%", call.pos, call.percent)
                        })
                        .collect::<Vec<_>>()
                        .join(", ")
                })
                .collect();
            format!("{} sn{}: {}", entry.kind, entry.sn, groups.join(" | "))
        })
        .collect();
    Some(format!("[{}]", sections.join("; ")))
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = "[DB sn1: \u{0394}5 100%, \u{0394}8 100%; OH sn1: 11 100%]";

    #[test]
    fn parses_the_documented_tail() {
        let parsed = parse_tail(SAMPLE).expect("well-formed tail");
        assert_eq!(parsed.len(), 2);
        assert_eq!(parsed[0].kind, "DB");
        assert_eq!(parsed[0].sn, 1);
        assert_eq!(
            parsed[0].alternatives,
            vec![vec![
                Call {
                    pos: 5,
                    percent: 100
                },
                Call {
                    pos: 8,
                    percent: 100
                }
            ]]
        );
        assert_eq!(parsed[1].kind, "OH");
        assert!(parsed[1].is_localized());
    }

    /// `|` separates mutually exclusive candidates for one feature — the
    /// engine could not choose between them.
    #[test]
    fn alternatives_are_kept_apart() {
        let parsed =
            parse_tail("[DB sn1: \u{0394}5 100% | \u{0394}14 50% | \u{0394}15 50%]").unwrap();
        assert_eq!(parsed[0].alternatives.len(), 3);
        assert!(parsed[0].is_localized(), "each group holds one candidate");

        let split = parse_tail("[DB sn1: \u{0394}14 50%, \u{0394}15 50%]").unwrap();
        assert_eq!(split[0].alternatives.len(), 1);
        assert!(!split[0].is_localized(), "one group, two candidates");
    }

    #[test]
    fn tokens_carry_no_whitespace() {
        let parsed = parse_tail(SAMPLE).unwrap();
        let tokens = tokens(&parsed, |_| None);
        assert_eq!(tokens, ["dbPos(sn1:5@100,8@100)", "mPos(sn1:11OH@100)"]);
        for token in &tokens {
            assert!(
                !token.contains(char::is_whitespace),
                "a `.smi` reader would truncate {token}"
            );
        }
    }

    /// A token about a floating group anchors to the stub's own label, so it
    /// survives the molecule being renumbered.
    #[test]
    fn a_stub_anchor_replaces_the_chain_anchor() {
        let parsed = parse_tail("[OH sn1: 11 50%, 13 50%]").unwrap();
        assert_eq!(
            tokens(&parsed, |_| Some("OH1".to_string())),
            ["mPos(OH1:11OH@50,13OH@50)"]
        );
    }

    #[test]
    fn tokens_round_trip_through_the_tail() {
        for tail in [
            SAMPLE,
            "[DB sn1: \u{0394}5 100%, \u{0394}8 100%, \u{0394}12 100% | \u{0394}14 50% | \u{0394}15 50%]",
            "[OH sn2: 11 50%, 13 50%]",
            "[DB sn1: \u{0394}9 92%; DB sn2: \u{0394}9 100%; OH sn2: 5 80%]",
        ] {
            let parsed = parse_tail(tail).unwrap_or_else(|| panic!("parse {tail}"));
            let encoded = tokens(&parsed, |_| None).join(";");
            let decoded = from_tokens(&encoded, |anchor| anchor.strip_prefix("sn")?.parse().ok());
            assert_eq!(decoded, parsed, "tokens lost something: {encoded}");
            assert_eq!(format_tail(&decoded).as_deref(), Some(tail), "{encoded}");
        }
    }

    #[test]
    fn malformed_tails_are_refused() {
        for tail in [
            "not a tail",
            "[DB sn1]",
            "[DB sn1: ]",
            "[DB: 5 100%]",
            "[DB sn1: five 100%]",
            "[DB snX: 5 100%]",
        ] {
            assert_eq!(parse_tail(tail), None, "{tail} should not parse");
        }
    }
}

#[cfg(test)]
mod round_trip {
    use crate::{name2smiles, smiles2name};

    /// A name's consensus tail must survive the trip through the structure —
    /// which is the whole point of carrying it in the trailer instead of
    /// throwing it away.
    #[test]
    fn consensus_survives_the_round_trip() {
        for name in [
            "FA 18:1(9) [DB sn1: \u{0394}9 92%]",
            "FA 18:2(9,12) [DB sn1: \u{0394}9 100%, \u{0394}12 88%]",
            "PC 16:0/18:1(9) [DB sn2: \u{0394}9 100%]",
            "PC 16:0_18:1(9) [DB sn2: \u{0394}9 100%]",
            // An unlocalized group: the token anchors to the stub's own label
            // and narrows the `m:` block from "any carbon" to two candidates.
            "FA 20:4;OH [OH sn1: 11 50%, 13 50%]",
        ] {
            let smiles = name2smiles(name).unwrap_or_else(|| panic!("{name} should resolve"));
            assert_eq!(smiles2name(&smiles).as_deref(), Some(name), "{name}");
        }
    }

    /// The consensus rides in the trailer and changes no CXSMILES field, so
    /// the structure is byte-identical with and without a tail.
    #[test]
    fn consensus_does_not_touch_the_structure() {
        let plain = name2smiles("FA 18:2(9,12)").unwrap();
        let annotated = name2smiles("FA 18:2(9,12) [DB sn1: \u{0394}9 100%]").unwrap();
        assert_eq!(
            plain.split(" |").next(),
            annotated.split(" |").next(),
            "the molecule changed"
        );
        assert!(!plain.contains("dbPos"), "{plain}");
        assert!(annotated.contains("dbPos(sn1:9@100)"), "{annotated}");
    }

    /// A tail this crate cannot parse is ignored rather than guessed at, and
    /// the structure comes out exactly as it would without one.
    #[test]
    fn an_unparsable_tail_is_ignored() {
        assert_eq!(
            name2smiles("FA 18:1(9) [not a consensus]"),
            name2smiles("FA 18:1(9)")
        );
    }
}
