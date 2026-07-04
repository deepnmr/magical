# MBP example (maltose-binding protein) — performance/scale demo only

A realistic-**scale** example: **192 methyl peaks / 142 methyl sites** built from
the real 3D structure of maltose-binding protein (PDB **1ANF**). Its purpose is to
exercise the pipeline at real protein size (~20 s runtime), **not** to produce
correct assignments.

> **⚠️ This is a synthetic dataset and its assignments are NOT meaningful.**
>
> - The 3D coordinates are real (1ANF), but the HMQC/NOESY peak lists are
>   **generated**: chemical shifts are type-based (a per-residue-type constant +
>   deterministic jitter), and NOE intensities follow a 1/r⁶ law.
> - Because the synthetic shifts of same-type methyls nearly overlap, the data
>   carries almost no information to distinguish, say, one Leu from another. The
>   pipeline recovers **only ~7% of the generation ground truth** (assignments are
>   type-correct but residue-arbitrary within a type).
> - It therefore **cannot** be compared to experimental assignments (e.g. BMRB),
>   and should not be read as biology. It is a *speed and scale* fixture only.
>
> A meaningful accuracy test needs **real** experimental methyl shifts (e.g. from
> a BMRB entry) plus a matching NOESY peak list — see "Real validation" below.

## Run it

```bash
python -m magic run examples/mbp/control.txt --output-dir mbp_run
head mbp_run/Output/assignments.tsv
cat mbp_run/Output/summary.json
```

Typical runtime: ~20 s. All 192 peaks get assigned.

## Files

| file | contents |
|---|---|
| `control.txt` | recipe (labeling `A;I,CD1;L,CD1,CD2;M;T;V,CG1,CG2`, cutoffs) |
| `hmqc.list` | 192 methyl peaks (¹³C/¹H shift + residue type) |
| `noesy.list` | 808 NOE cross-peaks from real inter-methyl distances (<6 Å) |
| `seq.txt` | 142 methyl-bearing residues |
| `model.pdb` | methyl-carbon coordinates from 1ANF |

To use the *full* 1ANF structure instead of the stripped methyl-carbon PDB,
download it from https://files.rcsb.org/download/1ANF.pdb and point `PDB =` at it.

## Real validation (not done here)

To check MAGIC against experimental truth you would:

1. Pull real methyl ¹³C/¹H shifts for MBP from a **BMRB** entry and build `hmqc.list`
   from them (with the correct residue types).
2. Provide a real methyl-methyl NOESY peak list (or simulate one from 1ANF using the
   **real** shifts, not type-based ones).
3. Run MAGIC and compare `best_assignment` to the BMRB residue label per peak.

This example does none of that; its shifts are fabricated, so its assignments are
by construction not comparable to BMRB.
