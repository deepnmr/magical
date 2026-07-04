"""Structure parsing and model-network construction."""

from collections import defaultdict
from math import dist
from typing import Dict, Iterable, List, Tuple

import numpy as np

from .models import MethylSite


THREE_TO_ONE = {
  'ALA': 'A',
  'ILE': 'I',
  'LEU': 'L',
  'MET': 'M',
  'THR': 'T',
  'VAL': 'V',
}

DEFAULT_LABELING = {
  'A': ('CB',),
  'I': ('CD1',),
  'L': ('CD1', 'CD2'),
  'M': ('CE',),
  'T': ('CG2',),
  'V': ('CG1', 'CG2'),
}


def parse_labeling(labeling: str) -> Dict[str, Tuple[str, ...]]:
  parsed: Dict[str, Tuple[str, ...]] = {}
  for chunk in labeling.split(';'):
    stripped = chunk.strip()
    if not stripped:
      continue
    tokens = [token.strip() for token in stripped.split(',') if token.strip()]
    residue = tokens[0]
    atoms = tuple(tokens[1:]) if len(tokens) > 1 else DEFAULT_LABELING.get(residue, ())
    if atoms:
      parsed[residue] = atoms
  return parsed


def parse_pdb_coordinates(pdb_lines: Iterable[str]) -> Dict[Tuple[str, int], Dict[str, Tuple[float, float, float]]]:
  coordinates: Dict[Tuple[str, int], Dict[str, Tuple[float, float, float]]] = defaultdict(dict)
  for line in pdb_lines:
    fields = line.split()
    if not fields or fields[0] != 'ATOM':
      continue
    residue_name = fields[3]
    if residue_name not in THREE_TO_ONE:
      continue
    residue_number = int(fields[5])
    atom_name = fields[2]
    coordinates[(THREE_TO_ONE[residue_name], residue_number)][atom_name] = (
      float(fields[6]),
      float(fields[7]),
      float(fields[8]),
    )
  return coordinates


def build_methyl_sites(
  sequence: Iterable[Tuple[str, int]],
  pdb_lines: Iterable[str],
  labeling: str,
) -> Tuple[MethylSite, ...]:
  labeling_map = parse_labeling(labeling)
  pdb_coords = parse_pdb_coordinates(pdb_lines)
  sites: List[MethylSite] = []
  for residue_type, residue_number in sequence:
    atoms = labeling_map.get(residue_type, ())
    if not atoms:
      continue
    residue_coords = pdb_coords.get((residue_type, residue_number), {})
    if residue_type in {'L', 'V'} and len(atoms) == 2:
      missing = [atom for atom in atoms if atom not in residue_coords]
      if missing:
        raise ValueError(f'Missing atoms for {residue_type}{residue_number}: {", ".join(missing)}')
      base_name = 'CD' if residue_type == 'L' else 'CG'
      sites.append(
        MethylSite(
          index=len(sites),
          label=f'{residue_type}{residue_number}{base_name}',
          residue_type=residue_type,
          residue_number=residue_number,
          atom_label=base_name,
          member_atoms=tuple(atoms),
          coordinates=tuple(residue_coords[atom] for atom in atoms),
          capacity=2,
        )
      )
      continue

    for atom in atoms:
      if atom not in residue_coords:
        raise ValueError(f'Missing atom for {residue_type}{residue_number}: {atom}')
      sites.append(
        MethylSite(
          index=len(sites),
          label=f'{residue_type}{residue_number}{atom}',
          residue_type=residue_type,
          residue_number=residue_number,
          atom_label=atom,
          member_atoms=(atom,),
          coordinates=(residue_coords[atom],),
          capacity=1,
        )
      )
  return tuple(sites)


def distance_weight(distance_value: float, low_cut: float, max_distance: float) -> float:
  if distance_value <= low_cut:
    return 1.0
  if distance_value > max_distance:
    return 0.0
  return max(0.0, (max_distance - distance_value) / (max_distance - low_cut))


def site_pair_weight(site_a: MethylSite, site_b: MethylSite, low_cut: float, max_distance: float) -> float:
  coordinates_a = list(site_a.coordinates)
  coordinates_b = list(site_b.coordinates)
  if site_a.index == site_b.index:
    if len(coordinates_a) < 2:
      return 0.0
    return sum(
      distance_weight(dist(coordinates_a[left], coordinates_a[right]), low_cut, max_distance)
      for left in range(len(coordinates_a))
      for right in range(left + 1, len(coordinates_a))
    )
  return sum(
    distance_weight(dist(coord_a, coord_b), low_cut, max_distance)
    for coord_a in coordinates_a
    for coord_b in coordinates_b
  )


def build_model_matrix(sites: Tuple[MethylSite, ...], low_cut: float, max_distance: float) -> np.ndarray:
  size = len(sites)
  matrix = np.zeros((size, size), dtype=float)
  for site_a in sites:
    for site_b in sites:
      matrix[site_a.index, site_b.index] = site_pair_weight(site_a, site_b, low_cut, max_distance)
  return matrix
