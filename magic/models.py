"""Shared data models for the clean MAGIC implementation."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


@dataclass(frozen=True)
class HMQCPeak:
  index: int
  peak_id: str
  carbon_shift: float
  proton_shift: float
  residue_types: Tuple[str, ...]
  locked: bool = False
  geminal_partners: Tuple[str, ...] = ()


@dataclass(frozen=True)
class NoePeak:
  peak_id: str
  donor_carbon: float
  ref_carbon: float
  ref_proton: float
  intensity: float
  raw_fields: Tuple[str, ...]


@dataclass(frozen=True)
class NoesyExperiment:
  name: str
  nuclei: Tuple[str, ...]
  tolerances: Tuple[float, ...]
  peaks: Tuple[NoePeak, ...]

  @property
  def donor_carbon_tolerance(self) -> float:
    return float(self.tolerances[0])

  @property
  def ref_carbon_tolerance(self) -> float:
    return float(self.tolerances[-2])

  @property
  def ref_proton_tolerance(self) -> float:
    return float(self.tolerances[-1])


@dataclass(frozen=True)
class MethylSite:
  index: int
  label: str
  residue_type: str
  residue_number: int
  atom_label: str
  member_atoms: Tuple[str, ...]
  coordinates: Tuple[Tuple[float, float, float], ...]
  capacity: int = 1


@dataclass
class ConnectionEvidence:
  acceptor_index: int
  donor_index: int
  confidence: float
  score: float
  shared_noes: int
  donor_count: int
  experiment_name: str


@dataclass
class NetworkData:
  confidence_matrix: np.ndarray
  score_matrix: np.ndarray
  density_matrix: np.ndarray
  adjacency_matrix: np.ndarray
  overlap_groups: Tuple[Tuple[int, ...], ...]
  strips_by_experiment: Dict[str, Tuple[Tuple[NoePeak, ...], ...]]
  connections: List[ConnectionEvidence] = field(default_factory=list)


@dataclass
class AssignmentHypothesis:
  mapping: Dict[int, int]
  score: float

  def canonical(self) -> Tuple[Tuple[int, int], ...]:
    return tuple(sorted(self.mapping.items()))


@dataclass
class ClusterResult:
  peaks: Tuple[int, ...]
  hypotheses: List[AssignmentHypothesis]


@dataclass(frozen=True)
class ParsedInputs:
  control_file: Path
  output_dir: Path
  input_dir: Path
  raw_output_dir: Path
  hmqc_peaks: Tuple[HMQCPeak, ...]
  noesy_experiments: Tuple[NoesyExperiment, ...]
  methyl_sites: Tuple[MethylSite, ...]
  low_cut: float
  max_distance: float
  score_factor: float
  labeling: str
