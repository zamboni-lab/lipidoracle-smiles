//! The slice of LipidOracle's name handling that the SMILES generator needs.
//!
//! Extracted from `src/nomenclature.rs` in the main LipidOracle crate. Only
//! the input-boundary helpers live here; scoring, ranking and the rest of the
//! nomenclature machinery stay upstream.

/// Splits a display name into its canonical part and its optional bracketed
/// confidence tail.
///
/// The EAD engines emit idlevel4 names carrying a per-position consensus in
/// square brackets, e.g.
/// `"FA 20:4;11OH [DB sn1: Δ5 100%, Δ8 100%; OH sn1: 11 100%]"`. The tail is
/// information the structure formats cannot represent at all (see
/// `dev/SMILES.md` §4.4), so the generator strips it and works from the
/// canonical name.
pub fn split_display_name(name: &str) -> (String, Option<String>) {
    let trimmed = name.trim();
    if let Some(start) = trimmed.find(" [") {
        if trimmed.ends_with(']') && start > 0 {
            return (
                trimmed[..start].trim_end().to_string(),
                Some(trimmed[start + 1..].to_string()),
            );
        }
    }
    (trimmed.to_string(), None)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn plain_name_has_no_tail() {
        assert_eq!(
            split_display_name("PC 16:0/18:1(9Z)"),
            ("PC 16:0/18:1(9Z)".to_string(), None)
        );
    }

    #[test]
    fn confidence_tail_is_split_off() {
        let (core, tail) = split_display_name("FA 20:4;11OH [DB sn1: \u{0394}5 100%]");
        assert_eq!(core, "FA 20:4;11OH");
        assert_eq!(tail.as_deref(), Some("[DB sn1: \u{0394}5 100%]"));
    }

    /// A bracket that isn't a trailing tail must be left alone — bracket
    /// atoms are ordinary SMILES, and a name could legitimately end in one.
    #[test]
    fn unterminated_bracket_is_not_a_tail() {
        let (core, tail) = split_display_name("FA 18:1 [oops");
        assert_eq!(core, "FA 18:1 [oops");
        assert_eq!(tail, None);
    }
}
