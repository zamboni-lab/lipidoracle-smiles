# Lipid Nomenclature

This document describes the nomenclature supported by the SMILES generator in
LipidOracle and how ambiguity is encoded in SMILES/CXSMILES output.

## Overview

The SMILES generator converts LipidOracle's lipid nomenclature into: -
**SMILES**: Simplified Molecular Input Line Entry System for well-defined
structures - **CXSMILES**: Extended SMILES with `Sg:` (variable-length run),
`m:` (position-variation bond) and `RG:` (R-group alternatives) blocks for
ambiguous structures

## Supported Lipid Classes

  | Class                    | Format   | Example                      | Notes                                              |
  | ------------------------ | -------- | ---------------------------- | -------------------------------------------------- |
  | Fatty Acid               | FA       | `FA 20:4(5,8,11,14)`         | Acyl chain with carboxyl                           |
  | Fatty Acid Amide         | NAE      | `NAE 20:4(5,8,11,14)`        | N-acyl ethanolamine                                |
  | Carnitine                | CAR      | `CAR 18:0`                   | Acyl-carnitine                                     |
  | Cholesterol Ester        | CE       | `CE 18:1(9Z)`                | Steroid with acyl chain                            |
  | Sterol                   | ST       | `ST`                         | Cholesterol backbone only                          |
  | Monoacylglycerol         | MG       | `MG 18:1(9)`                 | 1 acyl chain on glycerol                           |
  | Diacylglycerol           | DG       | `DG 16:0/18:1(9)`            | 2 acyl chains on glycerol                          |
  | Triacylglycerol          | TG       | `TG 16:0/18:1(9)/18:2(9,12)` | 3 acyl chains on glycerol                          |
  | Phosphatidylcholine      | PC       | `PC 16:0/18:1(9Z)`           | Glycerol + 2 chains + phosphocholine               |
  | Phosphatidylethanolamine | PE       | `PE 16:0/18:1(9)`            | Glycerol + 2 chains + phosphoethanolamine          |
  | Phosphatidylserine       | PS       | `PS 16:0/18:1(9)`            | Glycerol + 2 chains + phosphoserine                |
  | Phosphatidylglycerol     | PG       | `PG 16:0/18:1(9)`            | Glycerol + 2 chains + phosphoglycerol              |
  | Phosphatidylinositol     | PI       | `PI 16:0/18:1(9)`            | Glycerol + 2 chains + phosphoinositol              |
  | Phosphatidic Acid        | PA       | `PA 16:0/18:1(9)`            | Glycerol + 2 chains + phosphate                    |
  | Lyso-PC/PE/PS/etc        | LPC, LPE | `LPC 16:0`                   | Single acyl chain form (single-chain shorthand OK) |
  | Ceramide                 | Cer      | `Cer d18:1(4)/16:0`          | Sphingoid base + acyl chain                        |
  | Ceramide Phosphate       | CerP     | `CerP d18:1(4)/16:0`         | Ceramide + phosphate                               |
  | Sphingomyelin            | SM       | `SM d18:1(4)/16:0`           | Sphingoid base + acyl + phosphocholine             |
  | Hexosylceramide          | HexCer   | `HexCer d18:1(4)/16:0`       | Ceramide + hexose                                  |
  | Cardiolipin              | CL       | `CL 18:2_18:2_18:2_18:2`     | 4 acyl chains + central glycerol                   |

## Nomenclature Conventions

### Chain Notation: `C:DB(positions)`

Each acyl or sphingoid chain is specified as: - `C` = number of carbons - `DB` =
number of double bonds - `(positions)` = explicit double-bond positions
(optional, Delta-numbering from carboxyl)

**Examples:** - `16:0` → 16 carbons, 0 double bonds (saturated) - `18:1(9)` → 18
carbons, 1 double bond at position 9 - `20:4(5,8,11,14)` → 20 carbons, 4 double
bonds at positions 5, 8, 11, 14 - `18:2` → 18 carbons, 2 double bonds (positions
unknown, will generate placeholder)

### Separator Semantics: `/` vs `_`

The separator between chains determines how regiochemistry (sn-position
assignment) is encoded:

  | Separator         | Meaning                             | Example            | SMILES Output                            |
  | ----------------- | ----------------------------------- | ------------------ | ---------------------------------------- |
  | `/`               | sn-positions known                  | `PC 16:0/18:1(9Z)` | Single connected molecule, no `f:` block |
  | `_`               | sn-positions unknown (regioisomers) | `PC 16:0_18:1(9)`  | Separate wildcard fragments + `f:` block |
  | (none, shorthand) | Ambiguous composition               | `PC 19:2`          | **REJECTED** (returns None)              |

**Key rule:** Multi-chain lipids require explicit chains. Shorthand like
`PC 19:2` is rejected because it doesn't specify which chains have which
properties. Single-chain lipids (FA, CE, LPC) allow shorthand.

### Modifications: Functional Groups at Known Positions

Modifications are specified after a semicolon and can include multiple types:

  | Notation                    | Meaning                 | SMILES Result                      | Example                      |
  | --------------------------- | ----------------------- | ---------------------------------- | ---------------------------- |
  | `;posOH`                    | Hydroxyl groups         | `C(O)` at specified positions      | `FA 18:1(9);3OH,5OH`          |
  | `;posoxo`                   | Ketone groups           | `C(=O)` at specified positions     | `FA 18:0;5oxo`                |
  | `;COOH(pos)`                | Extra carboxyl branches | `C(C(=O)O)` at specified positions | `FA 20:4(5,8,11,14);COOH(8)` |
  | `;ep(pos)` / `;epox(pos)`   | Epoxide rings           | Parsed (rendering deferred)        | `FA 18:0;ep(5)`              |
  | `;cyc(pos)` / `;cyclo(pos)` | Cyclopropane rings      | Parsed (rendering deferred)        | `FA 18:0;cyc(3)`             |

### Double-Bond Geometry: `Z` vs `E` Notation

When double-bond position is known, you can optionally specify geometry:

  | Notation | Meaning                              | SMILES  | Example    |
  | -------- | ------------------------------------ | ------- | ---------- |
  | `(pos)`  | Position known, geometry unspecified | `C=C`   | `18:1(9)`  |
  | `(posZ)` | Position known, cis-double bond      | `/C=C\` | `18:1(9Z)` |
  | `(posE)` | Position known, trans-double bond    | `/C=C/` | `18:1(9E)` |

## How Ambiguity is Encoded in CXSMILES

### 1. Unlocalized Double Bonds → `Sg:` Block

When a chain declares double bonds but doesn't specify all positions:

**Input:** `FA 20:4` (20 carbons, 4 double bonds, positions unknown)

**Output:** SMILES with an unspecified double-bond run plus CXSMILES `Sg:` blocks
```
OC(=O)CC=CC=CC=CC=CC |Sg:n:3:a:ht,Sg:n:5:b:ht,Sg:n:7:c:ht,Sg:n:9:d:ht,Sg:n:12:e:ht| a+b+c+d+e=14
```

- Unlocalized positions are represented by an `Sg:` flexible run; the terminal
  carbon keeps the chain from ending in `C=C`
- Five markers for four double bonds (an `N+1` rule: one run before the first
  bond, then one per bond)
- The trailing `a+b+c+d+e=14` size constraint keeps the chain at 20 carbons. It
  is a LipidOracle/CDK-cookbook convention, not part of the CXSMILES standard
- Unspecified *geometry* needs no annotation at all: a plain `C=C` with no `/`
  or `\` already means "cis or trans, not determined". (Earlier versions added a
  `ctu:` block here; that block is a *query* feature for matching either
  configuration when searching, and added nothing to a structure.)

**Strategy:** This allows the structure to be rendered while clearly marking
which double bonds are not actually localized.

### 2. Unresolved Regiochemistry (sn-position) → `RG:` Block

When chains are joined with `_` (regioisomers) or when a lipid class always uses
`_`:

**Input:** `PC 16:0_18:1(9)` (unknown which chain is at sn-1 vs sn-2)

**Output:** Backbone with R-group slots + CXSMILES `RG:` alternatives
```
C(COP(=O)([O-])OCC[N+](C)(C)C)(O*)CO* |$;;;;;;;;;;;;;;R1;;;R1$,RG:_R1={C(=O)CCCCCCCCCCCCCCC},{C(=O)CCCCCCCC=CCCCCCCCC}|
```

**Structure breakdown:**

- **Backbone:** `C(COP(=O)([O-])OCC[N+](C)(C)C)(O*)CO*` with a `*` R-group slot
  at each of the two ester positions
- **Atom labels:** `$;;…;R1;;;R1$` — one `;`-separated slot per atom, in emission
  order, marking those two `*` atoms as `R1`
- **`RG:` block:** `RG:_R1={...},{...}` gives the two chains as the alternatives
  `R1` may take. A definition's first atom is implicitly its attachment point, so
  C1 of the chain bonds to the `O` carrying the `*`. The ester oxygen stays on the
  backbone, so each definition keeps all of its chain's own carbons.

**Strategy:** This represents the fact that we know the composition and
structures of individual chains, but not their assignment to sn-positions on the
glycerol backbone.

**Limit:** an `RG:` definition cannot contain a nested CXSMILES block (CDK
rejects it), so a chain that needs an `Sg:` or `m:` block of its own cannot
become an alternative. When any chain is in that position the chains are written
into the backbone in name order and no `RG:` block is emitted — the double-bond
ambiguity is kept and the sn ambiguity goes unexpressed. Earlier versions used an
`f:` component-group block here, which was wrong: `f:` joins components into one
entity (a salt, a hydrate) and expresses *and*, never *or*.

### 3. Multiple Ambiguities Combined

**Input:** `DG 18:1_18:2;5OH` (unknown sn-positions, both chains have
unlocalized double bonds, one has a localized hydroxyl)

**Output:**
```
C(CO)(OC(=O)CCCC(O)CC=CC=CC)COC(=O)CC=CC |Sg:n:11:a:ht,Sg:n:13:b:ht,Sg:n:16:c:ht,Sg:n:21:d:ht,Sg:n:24:e:ht| a+b+c=10,d+e=15
```

- Both chains need `Sg:`, so the sn ambiguity cannot also be expressed (see the
  limit in §2): the chains are esterified in place, in name order
- Hydroxyls at known positions are still rendered literally, as `C(O)`
- Each chain gets its own size constraint, comma-joined and matched to the
  markers by emission order rather than by variable name

## Examples by Ambiguity Type

### Fully Resolved (No Ambiguity)

```
FA 20:4(5,8,11,14)
→ OC(=O)CCCC=CCC=CCC=CCC=CCCCCC
   (no CXSMILES suffix needed)

PC 16:0/18:1(9Z)
→ C(COP(=O)([O-])OCC[N+](C)(C)C)(OC(=O)CCCCCCC/C=C\CCCCCCCC)COC(=O)CCCCCCCCCCCCCCC
   (sn-positions known → single connected molecule, / separator)
```

### Unlocalized Double Bonds Only

```
FA 18:2
→ OC(=O)CC=CC=CC |Sg:n:3:a:ht,Sg:n:5:b:ht,Sg:n:8:c:ht| a+b+c=14
   (two unlocalized double bonds, represented by three `Sg:` runs)

PC 16:0/18:2
→ C(COP(=O)([O-])OCC[N+](C)(C)C)(OC(=O)CC=CC=CC)COC(=O)CCCCCCCCCCCCCCC |Sg:n:16:a:ht,Sg:n:18:b:ht,Sg:n:21:c:ht| a+b+c=14
   (sn-positions known, but DB positions in 18:2 chain unknown)
```

### Unresolved Regiochemistry (sn-position)

```
DG 16:0_18:1(9)
→ C(CO)(O*)CO* |$;;;;R1;;;R1$,RG:_R1={C(=O)CCCCCCCCCCCCCCC},{C(=O)CCCCCCCC=CCCCCCCCC}|
   (both chains fully localized, so both become `R1` alternatives — ambiguous
    which is sn-1 vs sn-2)

TG 18:0_18:1_18:2
→ C(COC(=O)CC=CC=CC)(OC(=O)CC=CC)COC(=O)CCCCCCCCCCCCCCCCC |Sg:n:5:a:ht,Sg:n:7:b:ht,Sg:n:10:c:ht,Sg:n:14:d:ht,Sg:n:17:e:ht| a+b+c=14,d+e=15
   (three chains, sn-positions unknown, but 18:1 and 18:2 have unlocalized DBs —
    so `Sg:` wins and the sn ambiguity is not expressed)
```

### Unresolved Regiochemistry + Unlocalized Modifications

```
FA 20:4(5,8,11,14);3OH
→ OC(=O)CC(O)CC=CCC=CCC=CCC=CCCCCC
   (hydroxyl at position 3, all double-bond positions known → plain SMILES)

FA 18:0;OH
→ OC(=O)CCCCCCCCCCCCCCCCC.*O |m:20:3.4.5.6.7.8.9.10.11.12.13.14.15.16.17.18.19|
   (hydroxyl present but unpositioned → an `m:` position-variation bond. Atom 20
    is the `*`, not the oxygen: the variable end of the bond has to be a dummy
    atom carrying exactly one bond, or CDK ignores the block and RDKit rejects
    the string outright.)

DG 18:1(9);5OH_18:0;3OH
→ C(CO)(O*)CO* |$;;;;R1;;;R1$,RG:_R1={C(=O)CCCC(O)CCCC=CCCCCCCCC},{C(=O)CC(O)CCCCCCCCCCCCCCC}|
   (unresolved sn-positions, and each chain's modifications are localized, so
    both chains become `R1` alternatives)
```

## SMILES Validation Rules

The generated SMILES must satisfy:

1. **Balanced parentheses and brackets:** Every `(` has matching `)`, every `[`
   has matching `]`
2. **Even ring-closure digits:** Every numeric ring closure digit appears
   exactly twice
3. **Valid atom symbols:** Only organic atoms (C, N, O, P, S, etc.) and aromatic
   variants (c, n, o, p, s)
4. **Proper CXSMILES format:** Extensions use `|key:values|` syntax

All generated SMILES are validated before returning.

## Limitations and Unsupported Features

### Shorthand Multi-Chain Lipids Are Rejected

```
PC 19:2        → None (ambiguous: which chains have which DBs?)
TG 54:3        → None (ambiguous: how are 54 carbons/3 DBs distributed?)
```

**Reason:** Shorthand notation without explicit chain breakdown is inherently
ambiguous. Without knowing the individual chain compositions, there's no
canonical SMILES structure. Users must specify explicit chains.

**Solution:** Provide explicit chains:
```
PC 16:0/18:2   → accepted
PC 16:0_18:2   → accepted (regioisomers)
```

### Generic Oxygen Without Position Breakdown

```
FA 18:1(9);O     → None (ambiguous: hydroxyl? ketone? carboxyl?)
FA 18:1(9);O2    → None (ambiguous: where are 2 oxygens?)
```

**Reason:** Unlike double bonds (which have a standard placeholder convention at
C9, step-by-3), oxygen modifications have no universal fallback position scheme.
Guessing would risk generating incorrect structures.

**Solution:** Provide specific position breakdown:
```
FA 18:1(9);OH(3,5)     → accepted (two hydroxyls)
FA 18:1(9);oxo(5)      → accepted (one ketone)
FA 18:1(9);OH(3);oxo(5) → accepted (mixed)
```

### Ring Structures (Epoxides, Cyclopropanes)

These are **parsed but not yet rendered** in SMILES:
```
FA 18:0;ep(5)    → Parser accepts, returns None (rendering deferred)
FA 18:0;cyc(3)   → Parser accepts, returns None (rendering deferred)
```

**Reason:** Ring closure in SMILES requires explicit closure bonds with matching
digits, which adds structural complexity. Current focus is on acyclic
modifications (OH, oxo, COOH).

**Future enhancement:** Implement proper SMILES ring notation for epoxides and
cyclopropanes.

## LipidOracle Conventions for Expressing Double-Bond Position Consensus (idlevel4)

The idlevel4 annotation represents the highest level of structural assignment,
where individual double-bond positions are localized via EAD (Exhaustive
All-vs-All Decoupling) matching. The consensus of each position assignment is
explicitly encoded. It's an attempt to aggregate all top-matching C=C isomers
into a single description that encapsulates the consensus of a C=C being present
at any position of the acyl chains.

### Consensus Notation: `position~percentage`

Double-bond positions in idlevel4 are expressed with consensus scores using the
`~` separator:

**Format:** `C:DB(pos1~pct1,pos2~pct2,...)`

**Example interpretations:**

  | Notation               | Meaning                                             | Consensus           |
  | ---------------------- | --------------------------------------------------- | ------------------- |
  | `18:1(9~100%)`         | Single unambiguous position at C9                   | Certain (100%)      |
  | `18:1(9~80%,11~20%)`   | Position 9 with 80% consensus, position 11 with 20% | Ambiguous isomers   |
  | `20:4(5,8,11,14~100%)` | Positions 5,8,11,14 all known, 100% consensus       | All fully localized |
  | `20:4(5~100%,8,11,14)` | Position 5 is certain, others deduced               | Mixed certainty     |

### What the Percentage Represents

The percentage reflects the **relative abundance or assignment consensus** for
each position alternative:

- **100%**: This position is the single most likely assignment
- **50%**: Equally likely alternative to another position (e.g., regioisomers)
- **0-99%**: Partial consensus; other positions compete for the assignment

**Example scenario:**

```
DG 18:1(9~70%,11~30%)_18:0
```

This means:

- First chain: 18 carbons, 1 double bond
  - 70% likely at position 9 (most confident)
  - 30% likely at position 11 (alternative position)
- Second chain: 18 carbons, saturated
- Regiochemistry unresolved (sn-positions unknown, indicated by `_`)

### Aggregation Rules for idlevel4

When multiple idlevel3 candidates are aggregated into a single idlevel4 entry:

1. **Top scoring position wins:** The position with the highest EAD score gets
   100%
2. **Near-tie positions included:** Positions with scores within a defined
   cutoff (e.g., within top 1-2%) are included as alternatives
3. **Consensus reflects ranking:** Percentage = (position_score / top_score) ×
   100%

**Example aggregation:**

If EAD scoring produces these candidates: - `20:4(5,8,11,14)` --- EAD score 0.95 -
`20:4(5,8,11,15)` --- EAD score 0.92 - `20:4(5,8,12,14)` --- EAD score 0.85

The idlevel4 output might be:
```
20:4(5~100%,8~100%,11~100%,14~97%,15~97%,12~89%)
```

This shows that positions 5, 8, 11 are confidently assigned, while position 14
has a strong alternative at 15, and position 12 is a weaker third option.

### SMILES Rendering with Consensus

**Important:** SMILES and CXSMILES notation does **not** encode consensus
percentages directly. The SMILES represents the top-scoring structure only:

```
20:4(5~100%,8~100%,11~100%,14~97%,15~97%)
  ↓
OC(=O)CCCC=CCC=CCC=CCC=CCCCCCCCC
(renders position 14, the top candidate)
```

If you need to capture all alternatives with their consensus scores, use the
**full annotation CSV** (not the SMILES column), which preserves the complete
`~percentage` notation and includes all individual candidates as idlevel3.

### Typical Consensus Patterns

  | Pattern            | Meaning                             | Occurs When                                |
  | ------------------ | ----------------------------------- | ------------------------------------------ |
  | All 100%           | Fully consensus assignment          | All DB positions unambiguously matched     |
  | Mixed (>90%, <10%) | Clear winner with minor alternative | One position strongly preferred            |
  | 50%/50% split      | Regioisomers equally likely         | EAD scores are equal or very close         |
  | Multiple splits    | Ambiguous chain                     | Multiple equally-good position assignments |

### Consensus and Regiochemistry

Consensus percentages **apply per chain**, independent of sn-position
(regiochemistry) ambiguity:

```
PC 16:0/18:1(9~80%,11~20%)_16:1(7~100%)
  ↑                    ↑            ↑
sn-pos unknown (/_)    |            chain 3 with position consensus
                       └─ chain 2 with position consensus
```

This tells you: - sn-positions of all chains unknown (marked by `_`) - Within
the chains, position confidences are known per chain

## References

- **SMILES:** Weininger D. "SMILES, a chemical language and information system."
  J Chem Inf Comput Sci. 1988.
- **CXSMILES:** Dalby A, et al. "Description of several chemical structure file
  formats used by computer programs." J Chem Inf Comput Sci. 1992.
- **CXSMILES Format Specification:** ChemAxon. "ChemAxon Extended SMILES and
  SMARTS (CXSMILES and CXSMARTS)."
  https://docs.chemaxon.com/latest/formats_chemaxon-extended-smiles-and-smarts-cxsmiles-and-cxsmarts.html
- **Lipid Nomenclature:** IUPAC-IUBMB Lipid Nomenclature Standards
- **Modification Reference:** Fahy et al. "Lipid classification, structures and
  tools." PMC7707175 (Tables 1A-C)

## Troubleshooting

  | Issue                      | Cause                                 | Solution                                   |
  | -------------------------- | ------------------------------------- | ------------------------------------------ |
  | SMILES is empty            | Shorthand ambiguity or generic oxygen | Provide explicit chain breakdown           |
  | SMILES looks wrong         | Unlocalized DB positions              | Check `Sg:` blocks and the size constraint  |
  | Structure is disconnected  | Unresolved regiochemistry             | Check `f:` block for fragment grouping     |
  | Modification not appearing | Generic ;O without positions          | Specify ;OH(...), ;oxo(...), or ;COOH(...) |

--------------------------------------------------------------------------------

**Document Version:** 1.0

**Last Updated:** 2026-07-31

**Implementation:** src/smiles.rs in LipidOracle-RS
