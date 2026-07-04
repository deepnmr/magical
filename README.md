# magical

Clean-room implementation of the **MAGIC** methyl-NMR automated assignment
workflow (Monneau et al. 2017), with a modular package layout and an optimized,
output-preserving assignment search.

New here? Start with the **[step-by-step tutorial](TUTORIAL.md)**.

## Install

```bash
python3 -m pip install numpy psutil
```

## Usage

```bash
# run the pipeline from a control file
python -m magic run magic_control.txt --output-dir run_01

# validate a control file and the files it references
python -m magic validate magic_control.txt

# print or write a starter control file
python -m magic template magic_control.txt
```

`python Magic3v1.py ...` and the legacy positional form
`Magic3v1.py control.txt outdir` also work.

### Control file

Key/value format (relative paths resolve from the control file's location):

```
HMQC = hmqc.list
NOESY = noesy.list
PDB = model.pdb
SEQ = seq.txt
LABELING = A;I,CD1;L,CD1,CD2;M;T;V,CG1,CG2
GEMINAL = 1
CUTOFF_FACTOR = 1.0
DISTANCE_LIMITS = 7 10
SCORE_TOL_END = true
```

A realistic worked example (192 methyls, from PDB 1ANF) lives in
[`examples/mbp/`](examples/mbp/); a 3-line smoke-test example is in
[`tmp_synthetic/`](tmp_synthetic/).

## Package layout

```
magic/
  cli.py        command-line entry (argparse)
  config.py     control-file parsing + template
  inputs.py     path resolution / validation / output-dir naming
  io.py         HMQC / NOESY / sequence loaders
  models.py     dataclasses
  network.py    experimental peak-network construction
  structure.py  PDB parsing + model-network construction
  search.py     local & global assignment search (incremental scoring + B&B)
  output.py     TSV / JSON writers
  pipeline.py   end-to-end orchestration
```

## Changes & optimization

See [`CHANGES.md`](CHANGES.md) (also `CHANGES.pdf`) for the refactoring and the
output-preserving performance work, including computation-order comparison
flowcharts. Highlights on real protein structures:

| Dataset | methyl peaks | Original | After |
|---|---|---|---|
| 1UBQ | 43 | 143.5 s | 7.5 s |
| 1ANF (MBP) | 192 | did not finish | 19.2 s |
| 1D8C | 369 | did not finish | 21.9 s |

Output stays byte-identical (assignments); `best_score` differs by ≤ 2×10⁻¹³
from floating-point summation order.
