"""MAGIC control-file format: parsing and starter templates.

Two on-disk formats are supported:

* Key/value (``HMQC = ...``), the readable modern form.
* The original alternating-line format, where each value sits on a fixed
  line number.  Parsing falls back to this when the key/value parse fails.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import List


CONTROL_KEY_ALIASES = {
  'hmqc': 'hmqc_name',
  'hmqc_name': 'hmqc_name',
  'noesy': 'noesy_list',
  'noesy_list': 'noesy_list',
  'pdb': 'pdb_name',
  'pdb_name': 'pdb_name',
  'seq': 'seq_name',
  'sequence': 'seq_name',
  'seq_name': 'seq_name',
  'labeling': 'labeling',
  'geminal': 'flag_geminal',
  'flag_geminal': 'flag_geminal',
  'cutoff_factor': 'cutoff_factor',
  'distance_limits': 'distance_limits',
  'score_tol_end': 'score_tol_end',
  'isolated_pass': 'score_tol_end',
}

LEGACY_CONTROL_LINE_INDEX = {
  'hmqc_name': 2,
  'noesy_list': 4,
  'pdb_name': 6,
  'seq_name': 8,
  'labeling': 10,
  'flag_geminal': 13,
  'cutoff_factor': 15,
  'distance_limits': 17,
  'score_tol_end': 19,
}

REQUIRED_KEYS = (
  'hmqc_name',
  'noesy_list',
  'pdb_name',
  'seq_name',
  'labeling',
  'flag_geminal',
  'cutoff_factor',
  'distance_limits',
  'score_tol_end',
)

MAGIC_INPUT_TEMPLATE = """# MAGIC control file
# Key/value files are easier to read than the original alternating-line format.
# Relative paths are resolved from the location of this control file.
#
# Set SCORE_TOL_END = false to skip the final isolated-peak pass.

HMQC = example_hmqc.list
NOESY = example_cch_noesy.list
PDB = example_model.pdb
SEQ = example_sequence.txt
LABELING = A;I,CD1;L,CD1,CD2;M;T;V,CG1,CG2
GEMINAL = 1
CUTOFF_FACTOR = 1.0
DISTANCE_LIMITS = 6 10
SCORE_TOL_END = true
"""


@dataclass
class MagicControlFile:
  hmqc_name: str
  noesy_list: List[str]
  pdb_name: str
  seq_name: str
  labeling: str
  flag_geminal: float
  cutoff_factor: float
  low_cut: float
  max_distance: float
  score_tol_end: str


def _split_noesy(raw: str) -> List[str]:
  return [item.strip() for item in raw.split(';') if item.strip()]


def _parse_distance_limits(raw: str, error_message: str):
  parts = raw.split()
  if len(parts) != 2:
    raise ValueError(error_message)
  return float(parts[0]), float(parts[1])


def parse_key_value_control_file(lines) -> MagicControlFile:
  values = {}
  for raw_line in lines:
    stripped = raw_line.strip()
    if not stripped or stripped.startswith('#'):
      continue
    if '=' in stripped:
      key, value = stripped.split('=', 1)
    elif ':' in stripped:
      key, value = stripped.split(':', 1)
    else:
      continue
    canonical_key = CONTROL_KEY_ALIASES.get(key.strip().lower().replace('-', '_').replace(' ', '_'))
    if canonical_key is not None:
      values[canonical_key] = value.strip()

  missing = [key for key in REQUIRED_KEYS if key not in values]
  if missing:
    raise KeyError('Missing key/value fields: {}'.format(', '.join(missing)))

  low_cut, max_distance = _parse_distance_limits(
    values['distance_limits'],
    'DISTANCE_LIMITS must contain exactly two numbers: lowCut and max_distance.',
  )
  return MagicControlFile(
    hmqc_name=values['hmqc_name'],
    noesy_list=_split_noesy(values['noesy_list']),
    pdb_name=values['pdb_name'],
    seq_name=values['seq_name'],
    labeling=values['labeling'],
    flag_geminal=float(values['flag_geminal']),
    cutoff_factor=float(values['cutoff_factor']),
    low_cut=low_cut,
    max_distance=max_distance,
    score_tol_end=values['score_tol_end'],
  )


def parse_legacy_control_file(lines) -> MagicControlFile:
  if len(lines) <= max(LEGACY_CONTROL_LINE_INDEX.values()):
    raise ValueError(
      'Legacy control file is too short. Expected at least {} lines.'.format(
        max(LEGACY_CONTROL_LINE_INDEX.values()) + 1
      )
    )

  def line(key: str) -> str:
    return lines[LEGACY_CONTROL_LINE_INDEX[key]].strip()

  low_cut, max_distance = _parse_distance_limits(
    line('distance_limits'),
    'Legacy control file line 18 must contain lowCut and max_distance.',
  )
  return MagicControlFile(
    hmqc_name=line('hmqc_name'),
    noesy_list=_split_noesy(line('noesy_list')),
    pdb_name=line('pdb_name'),
    seq_name=line('seq_name'),
    labeling=line('labeling'),
    flag_geminal=float(line('flag_geminal')),
    cutoff_factor=float(line('cutoff_factor')),
    low_cut=low_cut,
    max_distance=max_distance,
    score_tol_end=line('score_tol_end'),
  )


def parse_magic_control_file(control_file) -> MagicControlFile:
  with open(control_file, 'r') as handle:
    lines = handle.readlines()
  try:
    return parse_key_value_control_file(lines)
  except (KeyError, ValueError):
    return parse_legacy_control_file(lines)


def emit_template(output_file=None):
  if output_file is None:
    print(MAGIC_INPUT_TEMPLATE.rstrip())
    return

  output_path = Path(output_file).expanduser()
  if output_path.exists():
    raise FileExistsError('Refusing to overwrite existing file: {}'.format(output_path))
  with open(output_path, 'w') as handle:
    handle.write(MAGIC_INPUT_TEMPLATE)
  print('Template written to {}'.format(output_path))
