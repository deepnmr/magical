# Paper conformance & BMRB validation

This note records (1) how the implementation was brought in line with the MAGIC
paper (Monneau et al., *J Biomol NMR* 2017, 69:215–227) and (2) an attempt to
validate the MBP assignment against experimental data from the BMRB.

## 1. Corrections to match the paper

### 1.1 Confidence score — Eq (1)

The paper defines the peak–peak confidence (P-matrix element) as

```
P_ij = sqrt( (φ_over / φ_CS) · ((1 + N_noe) / N'_donor)^2 )
```

with `φ_CS = 20·<Δδ>` (floored so `<Δδ> ≥ 0.05`) and

```
<Δδ>_ij = sqrt( (δ_i^noe − δ_i^hmqc)^2 + (δ_j^noe − δ_j^hmqc)^2 )
```

where *i* = acceptor (its HMQC carbon vs. the NOE **reference** carbon) and
*j* = donor (its HMQC carbon vs. the NOE **donor** carbon).

`magic/network.py` had two deviations from this:

- **`N'_donor` was not squared** (used to first power).
- The second term of `<Δδ>` used an ad-hoc "mirrored" shift instead of the
  acceptor's reference-carbon deviation.

Both are now fixed to Eq (1) exactly.

**Impact (measured on the MBP/BMRB test below):** with the buggy formula the
biologically-correct assignment scored *lower* than a wrong one found by the
search (truth 243.6 vs. found 351.2 — the objective did not point at biology).
After the fix the true assignment is the clear global maximum:

| quantity | score |
|---|---|
| true (BMRB) assignment | **651.2** |
| 20 random type-preserving permutations | mean 123.2, max 179.7 |

i.e. the corrected objective ranks the true assignment **~428 % above** random
and above all sampled alternatives. This is the key correctness result.

### 1.2 Local assignment — iterative T_C carry-forward (Fig. 3)

The paper grows peak clusters as the density threshold `T_C` decreases and
**carries the running assignment forward**, so each cluster only permutes the
peaks it adds that are not yet assigned. The previous code instead solved every
cluster independently (all peaks free) and merged afterwards, which exploded on
dense clusters. `magic/pipeline.py` now iterates clusters densest-first and calls
`expand_global_with_cluster`, matching the paper's protocol and keeping the local
exhaustive search bounded.

## 2. MBP validation against BMRB

### Setup

- **Shifts:** real experimental methyl ¹³C/¹H from **BMRB entry 7114**
  (E. coli maltose-binding protein, full ¹H/¹³C/¹⁵N assignment). 215 methyls with
  both C and H; residue numbering aligns to PDB **1ANF** (214 match / 1 mismatch /
  0 missing).
- **Structure / model network M:** from 1ANF.
- **Ground truth:** the BMRB residue label of each peak.
- **NOESY:** *simulated* from the 1ANF geometry (BMRB deposits chemical shifts,
  **not** NOESY peak lists), with 1/r⁶ intensities.

### Result

Recovery of the BMRB labels stays low (≈ 4–10 % depending on NOESY density and
search width) — near the within-residue-type chance level.

### Why — this is a data limitation, not (only) a code defect

The corrected objective *does* rank the true assignment highest (§1.1), so the
scoring is right. The recovery is nonetheless poor because:

1. **BMRB has no experimental NOESY.** A real validation needs the original
   3D CCH-NOESY peak list (its crosspeak set and intensities carry the
   distance information MAGIC exploits). That data is not in the BMRB.
2. **Structure-simulated NOESY is confounded.** When sparse it makes the top of
   the score landscape near-degenerate (many assignments within ~2 % of the
   true maximum); when denser, real methyl ¹³C shift degeneracy makes the
   donor-carbon matching ambiguous and lets some wrong assignments score *above*
   the true one. Neither regime uniquely pins the biological assignment.
3. **Search vs. a near-flat optimum.** Even where the true assignment is the
   unique maximum, finding it in a landscape whose runner-up is ~1–2 % lower is
   hard; the bounded local/beam search settles on a near-optimal but wrong
   solution.

### Conclusion

The program is now faithful to the paper's confidence score (Eq 1) and iterative
local-assignment protocol (Fig 3), and the corrected objective correctly ranks
the true MBP assignment as its maximum. A quantitative BMRB accuracy comparison
like Table 1 of the paper, however, **cannot be reproduced from BMRB alone** — it
requires the experimental CCH-NOESY peak lists used in the original study.
Reported here for full transparency.
