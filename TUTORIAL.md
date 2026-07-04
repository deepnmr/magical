# MAGIC — Step-by-Step Tutorial (Dummy Guide)

A hand-held walkthrough for a complete beginner. You will install the tool, run
the bundled example, understand every input and output file, then build your own
input from scratch. No prior knowledge of the codebase is assumed.

> **What does this tool do?** Given a protein 3D structure and methyl NMR peak
> lists (an HMQC "fingerprint" and a NOESY "who-is-near-whom" experiment), MAGIC
> figures out which peak belongs to which methyl group in the structure — the
> "assignment" problem. You give it peaks + a structure; it gives you a table
> saying *peak P1 = residue Ala1's methyl*, and so on.

---

## Step 0. Prerequisites

You need **Python 3.9+**. Check:

```bash
python3 --version
```

If that prints a version, you are good. If `python3` is not found, install Python
from https://www.python.org/downloads/ first.

---

## Step 1. Get the code

```bash
git clone https://github.com/deepnmr/magical.git
cd magical
```

You are now inside the project folder. Everything below is run from here.

---

## Step 2. Install the two dependencies

MAGIC needs `numpy` (math) and `psutil` (system info). Install them:

```bash
python3 -m pip install numpy psutil
```

If `pip` complains about permissions, either use a virtual environment
(recommended) or add `--user`:

```bash
# option A: virtual environment (isolated, cleanest)
python3 -m venv .venv
source .venv/bin/activate          # on Windows: .venv\Scripts\activate
python3 -m pip install numpy psutil

# option B: user install
python3 -m pip install --user numpy psutil
```

---

## Step 3. Run the bundled example (maltose-binding protein)

The repo ships a realistic example in `examples/mbp/` — **192 methyl peaks**
built from the real 3D structure of maltose-binding protein (PDB 1ANF). Run it:

```bash
python3 -m magic run examples/mbp/control.txt --output-dir first_run
```

It takes ~20 seconds. A new folder `first_run/` appears containing `Input/`
(a copy of what you fed in) and `Output/` (the results).

> **Heads-up:** `examples/mbp/` is a *synthetic* dataset — real coordinates from
> 1ANF, but the peak lists are generated from them (see
> [`examples/mbp/README.md`](examples/mbp/README.md)). It shows the tool at real
> scale; it is not experimental biology.
>
> Want an instant smoke test instead? There is also a tiny 3-residue example:
> `python3 -m magic run tmp_synthetic/control.txt --output-dir tiny_run`.

> **Note:** the output folder must not already exist — MAGIC refuses to
> overwrite. If you re-run, use a new `--output-dir` name (`first_run2`, …) or
> delete the old folder.

---

## Step 4. Read the results

Look at the main answer (first few rows):

```bash
head first_run/Output/assignments.tsv
```

```
peak_id  carbon_shift  proton_shift  best_assignment  alternatives  overlap_degree
P1       25.544        0.748         L135CD           L135CD        7
P2       24.588        0.676         L280CD           L280CD        10
P3       20.132        0.955         V50CG            V50CG         2
P4       19.176        0.883         V259CG           V259CG        5
P5       14.720        0.711         I108CD1          I108CD1       3
```

How to read it, column by column:

| column | meaning |
|---|---|
| `peak_id` | the peak name from your HMQC list |
| `carbon_shift`, `proton_shift` | that peak's ¹³C / ¹H chemical shifts |
| `best_assignment` | **the answer** — which methyl site the peak maps to (e.g. `L135CD` = Leu 135, atom CD) |
| `alternatives` | other plausible sites, if the data was ambiguous |
| `overlap_degree` | how many other peaks sit near this one (higher = more crowded region) |

So peak `P1` is Leu135's CD methyl, `P5` is Ile108's CD1, and so on.

Other output files (nice to know, not essential):

- `summary.json` — one-glance stats (how many peaks, how many assigned, best score).
- `connections.tsv` — the inferred "peak A is spatially near peak B" links, with confidence.
- `*_matrix.tsv` — the raw numeric matrices (score, confidence, density, model) for debugging.

Check the summary:

```bash
cat first_run/Output/summary.json
```

```json
{
  "num_peaks": 192,
  "num_methyl_sites": 142,
  "num_assigned_peaks": 192,
  "best_score": 1867.36,
  ...
}
```

`num_assigned_peaks` = 192 of 192 → every peak got assigned. 🎉

---

## Step 5. Understand the input files

Open `examples/mbp/` and look at each file (the tiny `tmp_synthetic/` files show
the exact same formats with just 3 lines each, if you prefer something you can
read at a glance). This is what you will replace with your own data later.

### 5a. `control.txt` — the recipe

```
HMQC = hmqc.list
NOESY = noesy.list
PDB = model.pdb
SEQ = seq.txt
LABELING = A;I,CD1;T
GEMINAL = 1
CUTOFF_FACTOR = 1.0
DISTANCE_LIMITS = 7 10
SCORE_TOL_END = true
```

- **HMQC / NOESY / PDB / SEQ** — filenames of the four data files (paths are
  relative to this control file).
- **LABELING** — which residue types are visible and which methyl atoms.
  Semicolon-separated. `A` means alanine with its default atom; `I,CD1` means
  isoleucine, atom CD1; `T` means threonine default. Full form for a real ILV
  sample: `A;I,CD1;L,CD1,CD2;M;T;V,CG1,CG2`.
- **GEMINAL** — 1 to treat the two methyls of Leu/Val as a linked pair, else 0.
- **CUTOFF_FACTOR** — scoring tolerance knob (1.0 is a good default).
- **DISTANCE_LIMITS** — two numbers `lowCut maxDistance` (Å). Methyl pairs closer
  than `lowCut` count fully; beyond `maxDistance` count as zero.
- **SCORE_TOL_END** — `true` to run the final isolated-peak pass, `false` to skip it.

### 5b. `hmqc.list` — the peak fingerprint

One line per methyl peak: `id  carbon_shift  proton_shift  residue_type(s)`.

```
P1 20.000 0.900 A
P2 13.000 0.800 I
P3 24.000 1.100 T
```

- Column 4 is the residue type(s) the peak could belong to. A single letter (`A`)
  or several stacked letters if uncertain (e.g. `IL` = could be Ile or Leu).
- Use `-` in column 4 to **lock** a peak to a known residue (advanced).
- An optional 5th column lists geminal partner peak ids (`;`-separated).

### 5c. `noesy.list` — the "who is near whom" experiment

```
13C;13C;1H            <- nuclei of each dimension
0.3;0.3;0.05          <- match tolerances per dimension
N1 13.000 20.000 0.900 100.0
N2 24.000 20.000 0.900 40.0
...
```

- **Line 1**: the nucleus of each column, `;`-separated.
- **Line 2**: the matching tolerance for each column (how close counts as "same").
- **Then one line per NOE cross-peak**: `id  donor_carbon  …  ref_carbon  ref_proton  intensity`.
  The first shift after the id is the *donor* carbon (the nearby methyl); the last
  three numbers are the *reference* carbon, reference proton, and the peak
  intensity. Bigger intensity = closer/stronger contact.

### 5d. `seq.txt` — the protein sequence

One residue per line, `<one-letter><number>`:

```
A1
I2
T3
```

Only methyl-bearing residues need to appear (A, I, L, M, T, V).

### 5e. `model.pdb` — the 3D structure

Standard PDB `ATOM` records. MAGIC only reads the methyl carbons of labeled
residues (plus needs the coordinates to compute distances):

```
ATOM      1  CB  ALA A   1       0.000   0.000   0.000  1.00 20.00           C
ATOM      3  CD1 ILE A   2       0.000   0.000   5.000  1.00 20.00           C
ATOM      5  CG2 THR A   3       0.000   0.000   9.000  1.00 20.00           C
```

You normally get this file straight from the RCSB PDB
(https://www.rcsb.org) for your protein.

---

## Step 6. Validate before you run (optional but smart)

Before a real run, check that your control file and all referenced files exist
and parse:

```bash
python3 -m magic validate examples/mbp/control.txt
```

```
Control file OK: .../examples/mbp/control.txt
HMQC: .../examples/mbp/hmqc.list
PDB:  .../examples/mbp/model.pdb
SEQ:  .../examples/mbp/seq.txt
NOESY: .../examples/mbp/noesy.list
```

If a file is missing or a field is malformed, it tells you here instead of
half-way through a run.

---

## Step 7. Make your own input

The fastest way to start a fresh control file:

```bash
python3 -m magic template my_control.txt
```

This writes a starter `my_control.txt` with every field and comments. Then:

1. Put your four data files next to it (`hmqc.list`, `noesy.list`, `model.pdb`,
   `seq.txt`) in the formats from Step 5.
2. Edit `my_control.txt` so the `HMQC =`, `NOESY =`, `PDB =`, `SEQ =` lines point
   at them, and set `LABELING` to match your isotope-labeling scheme.
3. Validate, then run:

```bash
python3 -m magic validate my_control.txt
python3 -m magic run my_control.txt --output-dir my_result
```

4. Read `my_result/Output/assignments.tsv` — same columns as Step 4.

---

## Step 8. Command cheat-sheet

```bash
python3 -m magic template [file]        # write / print a starter control file
python3 -m magic validate <control>     # check inputs without running
python3 -m magic run <control> --output-dir <dir>   # run the pipeline

# legacy entry points (identical behavior)
python3 Magic3v1.py run <control> --output-dir <dir>
python3 Magic3v1.py <control> <dir>     # shorthand: run + output dir
```

---

## Troubleshooting

| Symptom | Fix |
|---|---|
| `ModuleNotFoundError: No module named 'numpy'` | Run Step 2 (`pip install numpy psutil`). |
| `ModuleNotFoundError: No module named 'magic'` | Run from inside the project folder, or set `PYTHONPATH` to it. |
| `Output directory already exists` | Pick a new `--output-dir`, or delete the old one. |
| `Missing input files referenced by ...` | A path in `control.txt` is wrong; run `validate` to see which. |
| `Missing atom for X##: CD1` | Your PDB lacks that methyl atom, or `LABELING` names atoms your structure doesn't have. |
| A run takes very long on a big protein | Expected for dense/ambiguous data; the search is bounded by a node budget and will finish. See `CHANGES.md`. |

---

## Where to go next

- [`README.md`](README.md) — quick reference and package layout.
- [`CHANGES.md`](CHANGES.md) — how the assignment search is optimized, with flowcharts.
- Try a real structure: download a PDB from https://www.rcsb.org, write matching
  peak lists, and assign it.
