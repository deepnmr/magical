"""Input parsing for the clean MAGIC implementation."""

import shutil
from pathlib import Path
import re
from typing import List, Sequence, Tuple

from .models import HMQCPeak, NoePeak, NoesyExperiment


SEQUENCE_PATTERN = re.compile(r'([A-Z])(\d+)')


def load_hmqc(path: Path) -> Tuple[HMQCPeak, ...]:
  peaks: List[HMQCPeak] = []
  with open(path, 'r') as handle:
    for raw_line in handle:
      stripped = raw_line.strip()
      if not stripped or stripped.startswith('#'):
        continue
      fields = stripped.split()
      if len(fields) < 4:
        continue
      peak_id = fields[0]
      residue_types = tuple(sorted({token for token in fields[3] if token.isalpha()}))
      locked = fields[3] == '-'
      geminal_partners: Tuple[str, ...] = ()
      if len(fields) >= 5:
        raw_partners = [item for item in fields[4].split(';') if item and item != peak_id and item != '-']
        geminal_partners = tuple(raw_partners)
      peaks.append(
        HMQCPeak(
          index=len(peaks),
          peak_id=peak_id,
          carbon_shift=float(fields[1]),
          proton_shift=float(fields[2]),
          residue_types=residue_types,
          locked=locked,
          geminal_partners=geminal_partners,
        )
      )
  return tuple(peaks)


def load_noesy(path: Path) -> NoesyExperiment:
  with open(path, 'r') as handle:
    raw_lines = [line.strip() for line in handle if line.strip()]

  nuclei = tuple(token.strip() for token in raw_lines[0].split(';') if token.strip())
  tolerances = tuple(float(token.strip()) for token in raw_lines[1].split(';') if token.strip())
  peaks: List[NoePeak] = []
  for raw_line in raw_lines[2:]:
    fields = tuple(raw_line.split())
    if len(fields) < 5:
      continue
    peaks.append(
      NoePeak(
        peak_id=fields[0],
        donor_carbon=float(fields[1]),
        ref_carbon=float(fields[-3]),
        ref_proton=float(fields[-2]),
        intensity=float(fields[-1]),
        raw_fields=fields,
      )
    )
  return NoesyExperiment(name=path.name, nuclei=nuclei, tolerances=tolerances, peaks=tuple(peaks))


def load_sequence(path: Path) -> Tuple[Tuple[str, int], ...]:
  sequence: List[Tuple[str, int]] = []
  with open(path, 'r') as handle:
    for raw_line in handle:
      match = SEQUENCE_PATTERN.search(raw_line.strip())
      if match:
        sequence.append((match.group(1), int(match.group(2))))
  return tuple(sequence)


def copy_input_bundle(paths: Sequence[Path], target_directory: Path):
  target_directory.mkdir(parents=True, exist_ok=True)
  for source in paths:
    shutil.copy(str(source), str(target_directory / source.name))
