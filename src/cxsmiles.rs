//! Canonical SMILES serialization with lipid CXSMILES metadata renumbering.

use chematic_smiles::{canonical_atom_order, canonical_smiles, parse};

/// Canonicalizes a SMILES/CXSMILES string without detaching its lipid metadata.
///
/// CXSMILES atom numbers refer to the order in which atoms occur in the base
/// SMILES.  A normal canonical-SMILES rewrite therefore makes those numbers
/// stale.  This function uses the same canonical traversal for both jobs: it
/// writes the canonical base and applies the resulting old-to-new atom
/// permutation to:
///
/// * every `Sg:n:<atom>:...` flexible-chain marker;
/// * the floating atom and every candidate atom in `m:<atom>:<atoms>`;
/// * the positional `$...$` atom-label array used by `swappable(...)`;
/// * `atomProp:` atom indices.
///
/// Candidate indexes in each `m:` block are emitted in ascending canonical
/// atom order. Lipid trailer tokens after the closing pipe, such as `constrain(...)` and
/// `swappable(...)`, refer to names rather than atom positions and are retained
/// verbatim. Unknown CX fields are retained verbatim as well.
///
/// Returns `None` if the base SMILES is invalid, the CX block is unterminated,
/// or a recognized atom-indexed field is malformed or points outside the
/// molecule.
pub fn canonicalize(smi: &str) -> Option<String> {
    let ParsedCx {
        base,
        fields,
        trailer,
    } = split_cxsmiles(smi)?;

    // chematic-smiles currently rejects a disconnected component nested in a
    // branch and accepts only the bracketed wildcard spelling. Both occur in
    // lipid `m:` stubs. Normalize those spellings for parsing while retaining
    // the original-to-normalized atom permutation.
    let (lifted_base, old_to_parsed) = lift_nested_components(base)?;
    let parseable_base = bracket_bare_wildcards(&lifted_base);
    let mol = parse(&parseable_base).ok()?;
    let canonical = canonical_smiles(&mol);
    let canonical_mol = parse(&canonical).ok()?;

    // `canonical_atom_order` is canonical *rank* order, not necessarily the
    // textual order of branch atoms in the serialized SMILES. Pair the ranks
    // before and after serialization to obtain the actual atom permutation.
    let old_order = canonical_atom_order(&mol);
    let new_order = canonical_atom_order(&canonical_mol);
    if old_order.len() != new_order.len() {
        return None;
    }
    let mut parsed_to_new = vec![0; old_order.len()];
    for (old, new) in old_order.into_iter().zip(new_order) {
        parsed_to_new[old] = new;
    }
    let old_to_new = old_to_parsed
        .into_iter()
        .map(|parsed| parsed_to_new.get(parsed).copied())
        .collect::<Option<Vec<_>>>()?;

    let Some(fields) = fields else {
        return Some(canonical);
    };
    let rewritten = split_cx_fields(fields)
        .into_iter()
        .map(|field| rewrite_field(field, &old_to_new))
        .collect::<Option<Vec<_>>>()?;

    let mut out = canonical;
    out.push_str(" |");
    out.push_str(&rewritten.join(","));
    out.push('|');
    if !trailer.is_empty() {
        out.push(' ');
        out.push_str(trailer);
    }
    Some(out)
}

fn bracket_bare_wildcards(smiles: &str) -> String {
    let mut out = String::with_capacity(smiles.len());
    let mut in_brackets = false;
    for ch in smiles.chars() {
        match ch {
            '[' => {
                in_brackets = true;
                out.push(ch);
            }
            ']' => {
                in_brackets = false;
                out.push(ch);
            }
            '*' if !in_brackets => out.push_str("[*]"),
            _ => out.push(ch),
        }
    }
    out
}

/// Moves dot-disconnected components out of branches for the graph parser.
///
/// `A(B.*X)C` and `A(B)C.*X` describe the same two components, but moving the
/// stub changes textual atom order. The returned permutation maps each atom in
/// the original spelling to its index in the lifted spelling.
fn lift_nested_components(smiles: &str) -> Option<(String, Vec<usize>)> {
    let bytes = smiles.as_bytes();
    let mut ranges = Vec::<(usize, usize)>::new();
    let mut depth = 0usize;
    let mut in_brackets = false;
    let mut i = 0usize;
    while i < bytes.len() {
        match bytes[i] {
            b'[' => in_brackets = true,
            b']' => in_brackets = false,
            b'(' if !in_brackets => depth += 1,
            b')' if !in_brackets => depth = depth.checked_sub(1)?,
            b'.' if !in_brackets && depth > 0 => {
                let target_depth = depth;
                let mut scan_depth = depth;
                let mut scan_brackets = false;
                let mut end = i + 1;
                while end < bytes.len() {
                    match bytes[end] {
                        b'[' => scan_brackets = true,
                        b']' => scan_brackets = false,
                        b'(' if !scan_brackets => scan_depth += 1,
                        b')' if !scan_brackets => {
                            scan_depth = scan_depth.checked_sub(1)?;
                            if scan_depth < target_depth {
                                break;
                            }
                        }
                        _ => {}
                    }
                    end += 1;
                }
                if end == bytes.len() {
                    return None;
                }
                ranges.push((i, end));
                i = end;
                continue;
            }
            _ => {}
        }
        i += 1;
    }
    if in_brackets || depth != 0 {
        return None;
    }

    let atom_starts = atom_token_starts(smiles)?;
    let mut pieces = Vec::new();
    let mut cursor = 0usize;
    for &(start, end) in &ranges {
        pieces.push((cursor, start));
        cursor = end;
    }
    pieces.push((cursor, smiles.len()));
    pieces.extend(ranges.iter().copied());

    let mut lifted = String::with_capacity(smiles.len());
    let mut new_atoms_in_old_order = Vec::with_capacity(atom_starts.len());
    for (start, end) in pieces {
        lifted.push_str(&smiles[start..end]);
        new_atoms_in_old_order.extend(
            atom_starts
                .iter()
                .enumerate()
                .filter_map(|(old, &at)| (start <= at && at < end).then_some(old)),
        );
    }
    if new_atoms_in_old_order.len() != atom_starts.len() {
        return None;
    }
    let mut old_to_new = vec![0; atom_starts.len()];
    for (new, old) in new_atoms_in_old_order.into_iter().enumerate() {
        old_to_new[old] = new;
    }
    Some((lifted, old_to_new))
}

fn atom_token_starts(smiles: &str) -> Option<Vec<usize>> {
    let bytes = smiles.as_bytes();
    let mut starts = Vec::new();
    let mut i = 0usize;
    while i < bytes.len() {
        match bytes[i] {
            b'[' => {
                starts.push(i);
                i += 1;
                while i < bytes.len() && bytes[i] != b']' {
                    i += 1;
                }
                if i == bytes.len() {
                    return None;
                }
            }
            b'B' | b'C' | b'N' | b'O' | b'P' | b'S' | b'F' | b'I' | b'b' | b'c' | b'n' | b'o'
            | b'p' | b's' | b'*' => starts.push(i),
            _ => {}
        }
        i += 1;
    }
    Some(starts)
}

struct ParsedCx<'a> {
    base: &'a str,
    fields: Option<&'a str>,
    trailer: &'a str,
}

fn split_cxsmiles(smi: &str) -> Option<ParsedCx<'_>> {
    let trimmed = smi.trim();
    let Some(open) = trimmed.find('|') else {
        return Some(ParsedCx {
            base: trimmed,
            fields: None,
            trailer: "",
        });
    };
    let close = trimmed[open + 1..].find('|')? + open + 1;
    Some(ParsedCx {
        base: trimmed[..open].trim_end(),
        fields: Some(&trimmed[open + 1..close]),
        trailer: trimmed[close + 1..].trim(),
    })
}

/// Commas inside an atom-label block are data, not CX field separators.
fn split_cx_fields(fields: &str) -> Vec<&str> {
    let mut result = Vec::new();
    let mut start = 0;
    let mut in_labels = false;
    let mut escaped = false;
    for (i, ch) in fields.char_indices() {
        if escaped {
            escaped = false;
            continue;
        }
        match ch {
            '\\' => escaped = true,
            '$' => in_labels = !in_labels,
            ',' if !in_labels => {
                result.push(fields[start..i].trim());
                start = i + 1;
            }
            _ => {}
        }
    }
    result.push(fields[start..].trim());
    result
}

fn rewrite_field(field: &str, old_to_new: &[usize]) -> Option<String> {
    if field.starts_with('$') && field.ends_with('$') && field.len() >= 2 {
        return rewrite_labels(&field[1..field.len() - 1], old_to_new);
    }
    if let Some(rest) = field.strip_prefix("Sg:n:") {
        let (atom, suffix) = rest.split_once(':')?;
        return Some(format!("Sg:n:{}:{suffix}", remap_index(atom, old_to_new)?));
    }
    if let Some(rest) = field.strip_prefix("m:") {
        let (floating, candidates) = rest.split_once(':')?;
        let floating = remap_index(floating, old_to_new)?;
        let mut candidates = candidates
            .split('.')
            .map(|atom| remap_index(atom, old_to_new))
            .collect::<Option<Vec<_>>>()?;
        candidates.sort_unstable();
        let candidates = candidates
            .into_iter()
            .map(|index| index.to_string())
            .collect::<Vec<_>>();
        return Some(format!("m:{floating}:{}", candidates.join(".")));
    }
    if let Some(rest) = field.strip_prefix("atomProp:") {
        let entries = rest
            .split(':')
            .map(|entry| {
                let (atom, property) = entry.split_once('.')?;
                Some(format!("{}.{}", remap_index(atom, old_to_new)?, property))
            })
            .collect::<Option<Vec<_>>>()?;
        return Some(format!("atomProp:{}", entries.join(":")));
    }
    Some(field.to_string())
}

fn remap_index(text: &str, old_to_new: &[usize]) -> Option<usize> {
    old_to_new.get(text.parse::<usize>().ok()?).copied()
}

fn rewrite_labels(labels: &str, old_to_new: &[usize]) -> Option<String> {
    let labels: Vec<&str> = labels.split(';').collect();
    if labels.len() > old_to_new.len() && labels[old_to_new.len()..].iter().any(|s| !s.is_empty()) {
        return None;
    }
    let mut reordered = vec![""; old_to_new.len()];
    for (old, label) in labels.into_iter().take(old_to_new.len()).enumerate() {
        reordered[old_to_new[old]] = label;
    }
    while reordered.last() == Some(&"") {
        reordered.pop();
    }
    Some(format!("${}$", reordered.join(";")))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn canonicalizes_plain_smiles() {
        assert_eq!(canonicalize("OCC"), canonicalize("CCO"));
    }

    #[test]
    fn equivalent_atom_orders_produce_the_same_cxsmiles() {
        assert_eq!(
            canonicalize("OCC |Sg:n:0:a:ht,$oxygen;;tail$| constrain(a=1)"),
            canonicalize("CCO |Sg:n:2:a:ht,$tail;;oxygen$| constrain(a=1)")
        );
        assert_eq!(
            canonicalize("CC.*O |m:2:0.1|"),
            canonicalize("O*.CC |m:1:2.3|")
        );
    }

    #[test]
    fn remaps_sg_and_preserves_constraint() {
        let input = "OC(=O)CC=CC |Sg:n:3:a:ht,Sg:n:6:b:ht| constrain(a+b=15)";
        let out = canonicalize(input).unwrap();
        assert!(out.ends_with("| constrain(a+b=15)"), "{out}");

        let parsed = split_cxsmiles(&out).unwrap();
        let atom_count = parse(parsed.base).unwrap().atom_count();
        for field in split_cx_fields(parsed.fields.unwrap()) {
            if let Some(rest) = field.strip_prefix("Sg:n:") {
                let atom: usize = rest.split(':').next().unwrap().parse().unwrap();
                assert!(atom < atom_count, "{field} in {out}");
            }
        }
        assert_eq!(canonicalize(&out).as_deref(), Some(out.as_str()));
    }

    #[test]
    fn remaps_every_m_index() {
        let input = "OC(=O)CCCC.*O |m:7:3.4.5.6|";
        let out = canonicalize(input).unwrap();
        let parsed = split_cxsmiles(&out).unwrap();
        let field = parsed.fields.unwrap();
        let rest = field.strip_prefix("m:").unwrap();
        let (floating, candidates) = rest.split_once(':').unwrap();
        let atom_count = parse(parsed.base).unwrap().atom_count();
        assert!(floating.parse::<usize>().unwrap() < atom_count);
        let candidates = candidates
            .split('.')
            .map(|i| i.parse::<usize>().unwrap())
            .collect::<Vec<_>>();
        assert!(candidates.iter().all(|&i| i < atom_count));
        assert_eq!(canonicalize(&out).as_deref(), Some(out.as_str()));
    }

    #[test]
    fn sorts_m_candidates_after_canonicalizing_a_cyclized_chain() {
        let generated = crate::name2smiles("FA 20:2(5,8);[11-15cy5;13OH];OH").unwrap();
        let canonical = canonicalize(&generated).unwrap();
        let fields = split_cxsmiles(&canonical).unwrap().fields.unwrap();
        let m = split_cx_fields(fields)
            .into_iter()
            .find(|field| field.starts_with("m:"))
            .unwrap();
        assert_eq!(m, "m:1:2.3.4.5.6.7.8.9.10.11.12.13.14.15.16.17.21.22.23");
    }

    #[test]
    fn remaps_labels_and_keeps_swappable_names() {
        let input = "OCC |$left;;right$| swappable(left,right)";
        let out = canonicalize(input).unwrap();
        assert!(out.ends_with("| swappable(left,right)"), "{out}");
        assert!(out.contains("$"));
        assert_eq!(canonicalize(&out).as_deref(), Some(out.as_str()));
    }

    #[test]
    fn rejects_bad_references_instead_of_emitting_stale_indices() {
        assert_eq!(canonicalize("CC |m:9:0.1|"), None);
        assert_eq!(canonicalize("CC |Sg:n:x:a:ht|"), None);
        assert_eq!(canonicalize("CC |m:1:0.2|"), None);
    }
}
