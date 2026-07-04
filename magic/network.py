"""Experimental peak-network construction from HMQC and NOESY data."""

import math
from typing import Dict, List, Sequence, Tuple

import numpy as np

from .models import ConnectionEvidence, HMQCPeak, NetworkData, NoePeak, NoesyExperiment


def detect_overlaps(peaks: Sequence[HMQCPeak], carbon_tolerance: float, proton_tolerance: float) -> Tuple[Tuple[int, ...], ...]:
  overlap_groups: List[set[int]] = [set() for _ in peaks]
  for left in range(len(peaks)):
    for right in range(left + 1, len(peaks)):
      if (
        abs(peaks[left].carbon_shift - peaks[right].carbon_shift) <= carbon_tolerance
        and abs(peaks[left].proton_shift - peaks[right].proton_shift) <= proton_tolerance
      ):
        overlap_groups[left].add(right)
        overlap_groups[right].add(left)
  return tuple(tuple(sorted(group)) for group in overlap_groups)


def assign_noes_to_strips(peaks: Sequence[HMQCPeak], experiment: NoesyExperiment) -> Tuple[Tuple[NoePeak, ...], ...]:
  strips: List[List[NoePeak]] = [[] for _ in peaks]
  for peak in peaks:
    matched = [
      noe_peak
      for noe_peak in experiment.peaks
      if abs(noe_peak.ref_carbon - peak.carbon_shift) <= experiment.ref_carbon_tolerance
      and abs(noe_peak.ref_proton - peak.proton_shift) <= experiment.ref_proton_tolerance
    ]
    strips[peak.index] = matched
  return tuple(tuple(strip) for strip in strips)


def count_shared_noes(strip_left: Sequence[NoePeak], strip_right: Sequence[NoePeak], tolerance: float) -> int:
  donor_right = [peak.donor_carbon for peak in strip_right]
  shared = 0
  for peak in strip_left:
    if any(abs(peak.donor_carbon - donor_shift) <= tolerance for donor_shift in donor_right):
      shared += 1
  return shared


def mirrored_shift_delta(strip: Sequence[NoePeak], target_shift: float, tolerance: float) -> float:
  if not strip:
    return tolerance
  return min([abs(peak.donor_carbon - target_shift) for peak in strip] + [tolerance])


def build_single_experiment_matrices(
  peaks: Sequence[HMQCPeak],
  experiment: NoesyExperiment,
  overlap_groups: Tuple[Tuple[int, ...], ...],
) -> Tuple[np.ndarray, np.ndarray, Tuple[Tuple[NoePeak, ...], ...], List[ConnectionEvidence]]:
  size = len(peaks)
  confidence = np.zeros((size, size), dtype=float)
  score = np.zeros((size, size), dtype=float)
  strips = assign_noes_to_strips(peaks, experiment)
  shared_matrix = np.zeros((size, size), dtype=int)
  mirrored_matrix = np.full((size, size), experiment.donor_carbon_tolerance, dtype=float)

  for acceptor_index in range(size):
    for donor_index in range(size):
      if acceptor_index == donor_index:
        continue
      shared_matrix[acceptor_index, donor_index] = count_shared_noes(
        strips[acceptor_index],
        strips[donor_index],
        experiment.donor_carbon_tolerance,
      )
      mirrored_matrix[acceptor_index, donor_index] = mirrored_shift_delta(
        strips[donor_index],
        peaks[acceptor_index].carbon_shift,
        experiment.donor_carbon_tolerance,
      )

  evidence: List[ConnectionEvidence] = []
  for acceptor in peaks:
    strip = strips[acceptor.index]
    strip_total_intensity = sum(noe_peak.intensity for noe_peak in strip)
    if strip_total_intensity <= 0:
      continue
    for noe_peak in strip:
      donor_candidates = [
        donor.index
        for donor in peaks
        if donor.index != acceptor.index
        and abs(noe_peak.donor_carbon - donor.carbon_shift) <= experiment.donor_carbon_tolerance
      ]
      donor_count = max(1, len(donor_candidates))
      for donor_index in donor_candidates:
        delta_primary = abs(noe_peak.donor_carbon - peaks[donor_index].carbon_shift)
        delta_secondary = mirrored_matrix[acceptor.index, donor_index]
        delta_rms = max(0.05, math.sqrt(delta_primary ** 2 + delta_secondary ** 2))
        phi_cs = 20.0 * delta_rms
        shared_noes = int(shared_matrix[acceptor.index, donor_index])
        confidence_no_overlap = math.sqrt(((1 + shared_noes) ** 2) / (donor_count * phi_cs))
        confidence_with_overlap = 0.0 if overlap_groups[donor_index] else confidence_no_overlap
        score_value = confidence_no_overlap * (noe_peak.intensity / strip_total_intensity)

        confidence[acceptor.index, donor_index] = max(confidence[acceptor.index, donor_index], confidence_with_overlap)
        score[acceptor.index, donor_index] = max(score[acceptor.index, donor_index], score_value)
        evidence.append(
          ConnectionEvidence(
            acceptor_index=acceptor.index,
            donor_index=donor_index,
            confidence=confidence_with_overlap,
            score=score_value,
            shared_noes=shared_noes,
            donor_count=donor_count,
            experiment_name=experiment.name,
          )
        )
  return confidence, score, strips, evidence


def build_network_data(peaks: Sequence[HMQCPeak], experiments: Sequence[NoesyExperiment]) -> NetworkData:
  if not experiments:
    raise ValueError('At least one NOESY experiment is required.')

  overlap_groups = detect_overlaps(
    peaks,
    carbon_tolerance=max(experiment.ref_carbon_tolerance for experiment in experiments),
    proton_tolerance=max(experiment.ref_proton_tolerance for experiment in experiments),
  )

  confidence_matrices = []
  score_matrices = []
  strips_by_experiment: Dict[str, Tuple[Tuple[NoePeak, ...], ...]] = {}
  evidence: List[ConnectionEvidence] = []
  for experiment in experiments:
    experiment_confidence, experiment_score, strips, experiment_evidence = build_single_experiment_matrices(
      peaks,
      experiment,
      overlap_groups,
    )
    confidence_matrices.append(experiment_confidence)
    score_matrices.append(experiment_score)
    strips_by_experiment[experiment.name] = strips
    evidence.extend(experiment_evidence)

  confidence_matrix = np.mean(np.stack(confidence_matrices), axis=0)
  score_matrix = np.sum(np.stack(score_matrices), axis=0)
  adjacency = (confidence_matrix > 0).astype(float)
  np.fill_diagonal(adjacency, 1.0)
  density_matrix = adjacency @ adjacency
  np.fill_diagonal(density_matrix, 0.0)
  return NetworkData(
    confidence_matrix=confidence_matrix,
    score_matrix=score_matrix,
    density_matrix=density_matrix,
    adjacency_matrix=adjacency,
    overlap_groups=overlap_groups,
    strips_by_experiment=strips_by_experiment,
    connections=evidence,
  )
