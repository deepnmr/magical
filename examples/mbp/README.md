# MBP example (maltose-binding protein)

A realistic-scale example: **192 methyl peaks / 142 methyl sites** built from the
real 3D structure of maltose-binding protein (PDB **1ANF**).

> **Important — this is a *synthetic* dataset.** The 3D coordinates are real
> (taken from 1ANF), but the HMQC/NOESY peak lists were **generated** from those
> coordinates: chemical shifts are type-based with deterministic jitter, and NOE
> intensities follow a 1/r⁶ distance law. It reproduces the *scale and structure*
> of a real methyl-NMR dataset for demonstration and benchmarking, but the
> assignments it produces are **not experimentally validated biology**.

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
