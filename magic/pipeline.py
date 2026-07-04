"""End-to-end pipeline for the clean MAGIC implementation."""

from datetime import datetime
from pathlib import Path
from typing import Dict, List

from .inputs import normalize_output_directory_name, resolve_control_inputs
from .io import copy_input_bundle, load_hmqc, load_noesy, load_sequence
from .models import AssignmentHypothesis, ParsedInputs
from .network import build_network_data
from .output import write_outputs
from .search import (
  build_candidate_map,
  cluster_membership,
  cluster_size_threshold,
  connected_orphan_order,
  mapping_score,
  merge_hypothesis_sets,
  prune_hypotheses,
  solve_local_cluster,
)
from .structure import build_methyl_sites, build_model_matrix


def _copy_runtime_inputs(control_file: Path, resolved_inputs: Dict[str, object], input_dir: Path):
  files_to_copy = [control_file, resolved_inputs['hmqc'], resolved_inputs['pdb'], resolved_inputs['seq'], *resolved_inputs['noesy']]
  copy_input_bundle(files_to_copy, input_dir)


def _prepare_inputs(control_file: Path, output_dir: Path) -> ParsedInputs:
  config, resolved_inputs = resolve_control_inputs(control_file)
  input_dir = output_dir / 'Input'
  raw_output_dir = output_dir / 'Output'
  _copy_runtime_inputs(control_file, resolved_inputs, input_dir)

  hmqc_peaks = load_hmqc(resolved_inputs['hmqc'])
  sequence = load_sequence(resolved_inputs['seq'])
  with open(resolved_inputs['pdb'], 'r') as handle:
    pdb_lines = handle.readlines()
  methyl_sites = build_methyl_sites(sequence, pdb_lines, config.labeling)
  experiments = tuple(load_noesy(path) for path in resolved_inputs['noesy'])
  return ParsedInputs(
    control_file=control_file,
    output_dir=output_dir,
    input_dir=input_dir,
    raw_output_dir=raw_output_dir,
    hmqc_peaks=hmqc_peaks,
    noesy_experiments=experiments,
    methyl_sites=methyl_sites,
    low_cut=config.low_cut,
    max_distance=config.max_distance,
    score_factor=config.cutoff_factor,
    labeling=config.labeling,
  )


def _unique_thresholds(density_matrix) -> List[float]:
  thresholds = sorted({float(value) for value in density_matrix.flatten() if value >= 2.0}, reverse=True)
  return thresholds if thresholds else [2.0]


def _initial_hypotheses(
  candidate_map,
  peaks,
  score_matrix,
  model_matrix,
):
  locked_mapping = {}
  for peak in peaks:
    if not peak.locked:
      continue
    candidates = candidate_map[peak.index]
    if len(candidates) != 1:
      raise ValueError(f'Locked peak {peak.peak_id} resolved to {len(candidates)} model sites instead of one.')
    locked_mapping[peak.index] = candidates[0]
  return [AssignmentHypothesis(mapping=locked_mapping, score=mapping_score(locked_mapping, score_matrix, model_matrix))]


def _merge_clusters(
  global_hypotheses,
  local_cluster_results,
  site_capacities,
  score_matrix,
  model_matrix,
  total_peaks,
  score_factor,
):
  merged_signatures = set()
  current = list(global_hypotheses)
  for cluster_result in sorted(local_cluster_results, key=lambda result: (len(result.hypotheses), -len(result.peaks))):
    if not cluster_result.hypotheses:
      continue
    if cluster_result.peaks in merged_signatures:
      continue
    updated = merge_hypothesis_sets(
      current,
      cluster_result.hypotheses,
      site_capacities=site_capacities,
      score_matrix=score_matrix,
      model_matrix=model_matrix,
      total_peaks=total_peaks,
      score_factor=score_factor,
    )
    if updated:
      current = updated
      merged_signatures.add(cluster_result.peaks)
  return current


def _assign_orphans(
  global_hypotheses,
  peaks,
  candidate_map,
  site_capacities,
  score_matrix,
  model_matrix,
  total_peaks,
  score_factor,
):
  if not global_hypotheses:
    return []
  assigned = set(global_hypotheses[0].mapping.keys())
  orphan_indices = [peak.index for peak in peaks if peak.index not in assigned]
  for orphan_index in connected_orphan_order(score_matrix, assigned, orphan_indices):
    orphan_cluster = solve_local_cluster(
      (orphan_index,),
      candidate_map=candidate_map,
      site_capacities=site_capacities,
      score_matrix=score_matrix,
      model_matrix=model_matrix,
      total_peaks=total_peaks,
      score_factor=score_factor,
      max_keep=64,
    )
    global_hypotheses = merge_hypothesis_sets(
      global_hypotheses,
      orphan_cluster.hypotheses,
      site_capacities=site_capacities,
      score_matrix=score_matrix,
      model_matrix=model_matrix,
      total_peaks=total_peaks,
      score_factor=score_factor,
      max_keep=512,
    )
    if global_hypotheses:
      assigned = set(global_hypotheses[0].mapping.keys())
  return global_hypotheses


def run_magic(control_file, output_dir=None):
  control_path = Path(control_file).expanduser().resolve()
  runtime_output_dir = Path(normalize_output_directory_name(output_dir, datetime.now())).resolve()
  if runtime_output_dir.exists():
    raise FileExistsError(f'Output directory already exists: {runtime_output_dir}')
  runtime_output_dir.mkdir(parents=True, exist_ok=False)

  inputs = _prepare_inputs(control_path, runtime_output_dir)
  model_matrix = build_model_matrix(inputs.methyl_sites, inputs.low_cut, inputs.max_distance)
  network = build_network_data(inputs.hmqc_peaks, inputs.noesy_experiments)
  candidate_map = build_candidate_map(inputs.hmqc_peaks, inputs.methyl_sites)
  site_capacities = {site.index: site.capacity for site in inputs.methyl_sites}
  global_hypotheses = _initial_hypotheses(candidate_map, inputs.hmqc_peaks, network.score_matrix, model_matrix)

  local_results = {}
  for tc_value in _unique_thresholds(network.density_matrix):
    threshold_size = cluster_size_threshold(tc_value)
    for seed_index in range(len(inputs.hmqc_peaks)):
      cluster = cluster_membership(network.adjacency_matrix, network.density_matrix, seed_index, tc_value)
      if len(cluster) < threshold_size:
        continue
      if cluster not in local_results:
        local_results[cluster] = solve_local_cluster(
          cluster,
          candidate_map=candidate_map,
          site_capacities=site_capacities,
          score_matrix=network.score_matrix,
          model_matrix=model_matrix,
          total_peaks=len(inputs.hmqc_peaks),
          score_factor=inputs.score_factor,
        )

  global_hypotheses = _merge_clusters(
    global_hypotheses,
    list(local_results.values()),
    site_capacities=site_capacities,
    score_matrix=network.score_matrix,
    model_matrix=model_matrix,
    total_peaks=len(inputs.hmqc_peaks),
    score_factor=inputs.score_factor,
  )
  global_hypotheses = _assign_orphans(
    global_hypotheses,
    inputs.hmqc_peaks,
    candidate_map,
    site_capacities=site_capacities,
    score_matrix=network.score_matrix,
    model_matrix=model_matrix,
    total_peaks=len(inputs.hmqc_peaks),
    score_factor=inputs.score_factor,
  )
  final_hypotheses = prune_hypotheses(
    global_hypotheses,
    total_peaks=len(inputs.hmqc_peaks),
    score_factor=inputs.score_factor,
    max_keep=256,
  )
  write_outputs(inputs, inputs.methyl_sites, model_matrix, network, final_hypotheses)
  return inputs.output_dir
