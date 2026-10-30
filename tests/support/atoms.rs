// Shared by the unit tests (`src/forward.rs`) and the integration tests
// (`tests/testdata.rs`), which both include it by path.

/// The atom symbols of a SMILES string, in order: bracket atoms whole,
/// organic-subset atoms one letter each.
pub(crate) fn atom_symbols(smiles: &str) -> Vec<String> {
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
                out.push(chars[start..=i.min(chars.len() - 1)].iter().collect());
                i += 1;
            }
            c @ ('B' | 'C' | 'N' | 'O' | 'P' | 'S' | 'F' | 'I' | 'c' | 'n' | 'o' | 's' | 'p'
            | '*') => {
                out.push(c.to_string());
                i += 1;
            }
            _ => i += 1,
        }
    }
    out
}
