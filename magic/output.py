"""Output writers for the clean MAGIC implementation."""

import json
from pathlib import Path
from typing import Dict, Iterable, Sequence

import numpy as np

from .models import AssignmentHypothesis, HMQCPeak, MethylSite, NetworkData, ParsedInputs


def write_numeric_matrix(path: Path, matrix: np.ndarray):
  np.savetxt(path, matrix, fmt='%.6f', delimiter='\t')


def write_assignments(
  output_dir: Path,
  peaks: Sequence[HMQCPeak],
  sites: Sequence[MethylSite],
  hypotheses: Sequence[AssignmentHypothesis],
  network: NetworkData,
):
  site_labels = {site.index: site.label for site in sites}
  best_mapping = hypotheses[0].mapping if hypotheses else {}
  alternatives: Dict[int, list[str]] = {}
  for hypothesis in hypotheses[: min(50, len(hypotheses))]:
    for peak_index, site_index in hypothesis.mapping.items():
      alternatives.setdefault(peak_index, [])
      label = site_labels[site_index]
      if label not in alternatives[peak_index]:
        alternatives[peak_index].append(label)

  assignments_path = output_dir / 'assignments.tsv'
  with open(assignments_path, 'w') as handle:
    handle.write('peak_id\tcarbon_shift\tproton_shift\tbest_assignment\talternatives\toverlap_degree\n')
    for peak in peaks:
      best_assignment = site_labels.get(best_mapping.get(peak.index, -1), 'UNASSIGNED')
      handle.write(
        f'{peak.peak_id}\t{peak.carbon_shift:.3f}\t{peak.proton_shift:.3f}\t{best_assignment}\t'
        f'{";".join(alternatives.get(peak.index, []))}\t{len(network.overlap_groups[peak.index])}\n'
      )


def write_connections(output_dir: Path, peaks: Sequence[HMQCPeak], sites: Sequence[MethylSite], network: NetworkData):
  path = output_dir / 'connections.tsv'
  with open(path, 'w') as handle:
    handle.write('acceptor\tdonor\tconfidence\tscore\tshared_noes\tdonor_count\texperiment\n')
    for connection in network.connections:
      handle.write(
        f'{peaks[connection.acceptor_index].peak_id}\t{peaks[connection.donor_index].peak_id}\t'
        f'{connection.confidence:.6f}\t{connection.score:.6f}\t{connection.shared_noes}\t'
        f'{connection.donor_count}\t{connection.experiment_name}\n'
      )


def write_summary(
  output_dir: Path,
  inputs: ParsedInputs,
  sites: Sequence[MethylSite],
  network: NetworkData,
  hypotheses: Sequence[AssignmentHypothesis],
):
  best_mapping = hypotheses[0].mapping if hypotheses else {}
  summary = {
    'control_file': str(inputs.control_file),
    'num_peaks': len(inputs.hmqc_peaks),
    'num_methyl_sites': len(sites),
    'num_noesy_experiments': len(inputs.noesy_experiments),
    'num_final_hypotheses': len(hypotheses),
    'num_assigned_peaks': len(best_mapping),
    'best_score': hypotheses[0].score if hypotheses else 0.0,
    'mean_overlap_degree': float(np.mean([len(group) for group in network.overlap_groups])) if network.overlap_groups else 0.0,
  }
  with open(output_dir / 'summary.json', 'w') as handle:
    json.dump(summary, handle, indent=2)


def write_outputs(
  inputs: ParsedInputs,
  sites: Sequence[MethylSite],
  model_matrix: np.ndarray,
  network: NetworkData,
  hypotheses: Sequence[AssignmentHypothesis],
):
  output_dir = inputs.raw_output_dir
  output_dir.mkdir(parents=True, exist_ok=True)
  write_numeric_matrix(output_dir / 'model_matrix.tsv', model_matrix)
  write_numeric_matrix(output_dir / 'confidence_matrix.tsv', network.confidence_matrix)
  write_numeric_matrix(output_dir / 'score_matrix.tsv', network.score_matrix)
  write_numeric_matrix(output_dir / 'density_matrix.tsv', network.density_matrix)
  write_assignments(output_dir, inputs.hmqc_peaks, sites, hypotheses, network)
  write_connections(output_dir, inputs.hmqc_peaks, sites, network)
  write_summary(output_dir, inputs, sites, network, hypotheses)
