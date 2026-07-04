"""Local and global assignment search for the clean MAGIC implementation."""

from collections import Counter
import re
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np

from .models import AssignmentHypothesis, ClusterResult, HMQCPeak, MethylSite


_RESIDUE_NUMBER_PATTERN = re.compile(r'\w(\d+)')


def generic_residue_id(label: str) -> str:
  """Extract the residue number from a peak/site label, e.g. ``I12`` -> ``12``."""

  return _RESIDUE_NUMBER_PATTERN.search(label).group(1)


def cluster_size_threshold(tc_value: float) -> int:
  if tc_value >= 5:
    return 5
  if tc_value >= 3:
    return 4
  return 3


def build_candidate_map(peaks: Sequence[HMQCPeak], sites: Sequence[MethylSite]) -> Dict[int, Tuple[int, ...]]:
  sites_by_type: Dict[str, List[int]] = {}
  for site in sites:
    sites_by_type.setdefault(site.residue_type, []).append(site.index)

  candidates: Dict[int, Tuple[int, ...]] = {}
  for peak in peaks:
    if peak.locked:
      locked_matches = [site.index for site in sites if site.label == peak.peak_id]
      if not locked_matches:
        residue_number = int(generic_residue_id(peak.peak_id))
        locked_matches = [
          site.index
          for site in sites
          if site.residue_type == peak.peak_id[0] and site.residue_number == residue_number
        ]
      if not locked_matches:
        raise ValueError(f'Locked peak {peak.peak_id} does not match any model methyl site.')
      candidates[peak.index] = tuple(locked_matches)
      continue

    peak_candidates: List[int] = []
    for residue_type in peak.residue_types:
      peak_candidates.extend(sites_by_type.get(residue_type, ()))
    if not peak_candidates:
      raise ValueError(f'Peak {peak.peak_id} has no candidate model sites for types {peak.residue_types}.')
    candidates[peak.index] = tuple(sorted(set(peak_candidates)))
  return candidates


def mapping_score(mapping: Dict[int, int], score_matrix: np.ndarray, model_matrix: np.ndarray) -> float:
  total = 0.0
  for acceptor_index, acceptor_site in mapping.items():
    for donor_index, donor_site in mapping.items():
      total += float(score_matrix[acceptor_index, donor_index] * model_matrix[acceptor_site, donor_site])
  return total


def relative_score_slack(score_factor: float, assigned_ratio: float) -> float:
  return min(0.5, 0.1 + 0.35 * score_factor * max(0.0, 1.0 - assigned_ratio))


def deduplicate_hypotheses(hypotheses: Iterable[AssignmentHypothesis]) -> List[AssignmentHypothesis]:
  best_by_mapping: Dict[Tuple[Tuple[int, int], ...], AssignmentHypothesis] = {}
  for hypothesis in hypotheses:
    canonical = hypothesis.canonical()
    if canonical not in best_by_mapping or hypothesis.score > best_by_mapping[canonical].score:
      best_by_mapping[canonical] = hypothesis
  return list(best_by_mapping.values())


def prune_hypotheses(
  hypotheses: Iterable[AssignmentHypothesis],
  total_peaks: int,
  score_factor: float,
  max_keep: int,
) -> List[AssignmentHypothesis]:
  unique_hypotheses = deduplicate_hypotheses(hypotheses)
  if not unique_hypotheses:
    return []
  best_score = max(hypothesis.score for hypothesis in unique_hypotheses)
  assigned_ratio = max(len(hypothesis.mapping) / max(1, total_peaks) for hypothesis in unique_hypotheses)
  slack = relative_score_slack(score_factor, assigned_ratio)
  if best_score > 0:
    cutoff = best_score * (1.0 - slack)
  else:
    cutoff = best_score - max(1.0, slack)
  survivors = [hypothesis for hypothesis in unique_hypotheses if hypothesis.score >= cutoff]
  survivors.sort(key=lambda item: (-item.score, len(item.mapping)))
  return survivors[:max_keep]


def cluster_membership(adjacency: np.ndarray, density: np.ndarray, seed_index: int, tc_value: float) -> Tuple[int, ...]:
  members = [seed_index]
  for candidate in range(adjacency.shape[0]):
    if candidate == seed_index:
      continue
    if adjacency[seed_index, candidate] > 0 and density[seed_index, candidate] > tc_value:
      members.append(candidate)
  return tuple(sorted(set(members)))


def enumerate_cluster_assignments(
  cluster_peaks: Tuple[int, ...],
  candidate_map: Dict[int, Tuple[int, ...]],
  site_capacities: Dict[int, int],
) -> List[Dict[int, int]]:
  ordered_peaks = sorted(cluster_peaks, key=lambda peak_index: (len(candidate_map[peak_index]), peak_index))
  assignments: List[Dict[int, int]] = []

  def backtrack(depth: int, current_mapping: Dict[int, int], usage: Counter):
    if depth == len(ordered_peaks):
      assignments.append(dict(current_mapping))
      return
    peak_index = ordered_peaks[depth]
    for site_index in candidate_map[peak_index]:
      if usage[site_index] >= site_capacities[site_index]:
        continue
      current_mapping[peak_index] = site_index
      usage[site_index] += 1
      backtrack(depth + 1, current_mapping, usage)
      usage[site_index] -= 1
      del current_mapping[peak_index]

  backtrack(0, {}, Counter())
  return assignments


def _keep_cutoff(best_score: float, slack: float) -> float:
  """Score threshold that ``prune_hypotheses`` will apply for a given best."""

  if best_score > 0:
    return best_score * (1.0 - slack)
  return best_score - max(1.0, slack)


SOLVE_NODE_BUDGET = 1_000_000


def solve_local_cluster(
  cluster_peaks: Tuple[int, ...],
  candidate_map: Dict[int, Tuple[int, ...]],
  site_capacities: Dict[int, int],
  score_matrix: np.ndarray,
  model_matrix: np.ndarray,
  total_peaks: int,
  score_factor: float,
  max_keep: int = 256,
  node_budget: int = SOLVE_NODE_BUDGET,
) -> ClusterResult:
  """Enumerate cluster assignments with branch-and-bound pruning.

  Scores are accumulated incrementally as peaks are placed (avoiding an
  O(k^2) rescore per complete assignment), and partial assignments whose
  optimistic upper bound cannot reach the keep-band that ``prune_hypotheses``
  applies are pruned.  Because the bound is admissible (all score/model terms
  are non-negative) and the running-best cutoff never exceeds the final
  cutoff, this yields the same surviving hypotheses as full enumeration.

  ``node_budget`` caps the number of complete assignments explored so a
  pathologically dense, poorly-differentiated cluster degrades gracefully
  (returning the best band found so far) instead of hanging.  Realistic
  clusters finish far below the budget, leaving their result unchanged.
  """

  ordered_peaks = sorted(cluster_peaks, key=lambda peak_index: (len(candidate_map[peak_index]), peak_index))
  if not ordered_peaks:
    hypotheses = [AssignmentHypothesis(mapping={}, score=0.0)]
    return ClusterResult(
      peaks=cluster_peaks,
      hypotheses=prune_hypotheses(hypotheses, total_peaks=total_peaks, score_factor=score_factor, max_keep=max_keep),
    )

  # Compact the relevant score/model sub-matrices into plain Python lists so
  # the hot inner loop avoids per-element numpy scalar indexing.
  site_ids = sorted({site for peak in ordered_peaks for site in candidate_map[peak]})
  site_local = {site: local for local, site in enumerate(site_ids)}
  peak_score = score_matrix[np.ix_(ordered_peaks, ordered_peaks)].tolist()
  site_model = model_matrix[np.ix_(site_ids, site_ids)].tolist()
  candidates_local = [tuple(site_local[site] for site in candidate_map[peak]) for peak in ordered_peaks]
  capacities_local = [site_capacities[site] for site in site_ids]

  # Admissible per-peak upper bound on the score any remaining peak can add.
  max_model = max((max(row) for row in site_model), default=0.0)
  cross_sum = [
    sum(peak_score[i][j] + peak_score[j][i] for j in range(len(ordered_peaks)) if j != i)
    for i in range(len(ordered_peaks))
  ]
  peak_bound = [
    max((peak_score[i][i] * site_model[ls][ls] for ls in candidates_local[i]), default=0.0)
    + cross_sum[i] * max_model
    for i in range(len(ordered_peaks))
  ]

  assigned_ratio = len(ordered_peaks) / max(1, total_peaks)
  slack = relative_score_slack(score_factor, assigned_ratio)

  depth_count = len(ordered_peaks)
  placed_local: List[int] = []
  hypotheses: List[AssignmentHypothesis] = []
  usage = [0] * len(site_ids)
  best_score = [0.0]
  leaves = [0]
  stop = [False]

  def backtrack(depth: int, running_score: float, remaining_bound: float, mapping: Dict[int, int]):
    if depth == depth_count:
      best_score[0] = max(best_score[0], running_score)
      hypotheses.append(AssignmentHypothesis(mapping=dict(mapping), score=running_score))
      leaves[0] += 1
      if leaves[0] >= node_budget:
        stop[0] = True
      return
    peak_index = ordered_peaks[depth]
    next_bound = remaining_bound - peak_bound[depth]
    self_term = peak_score[depth][depth]
    for local_site in candidates_local[depth]:
      if usage[local_site] >= capacities_local[local_site]:
        continue
      model_row = site_model[local_site]
      delta = self_term * model_row[local_site]
      for placed_depth in placed_local:
        placed_site = mapping_local[placed_depth]
        delta += peak_score[depth][placed_depth] * model_row[placed_site]
        delta += peak_score[placed_depth][depth] * site_model[placed_site][local_site]
      new_score = running_score + delta
      if new_score + next_bound < _keep_cutoff(best_score[0], slack):
        continue
      usage[local_site] += 1
      mapping[peak_index] = site_ids[local_site]
      mapping_local[depth] = local_site
      placed_local.append(depth)
      backtrack(depth + 1, new_score, next_bound, mapping)
      placed_local.pop()
      usage[local_site] -= 1
      del mapping[peak_index]
      if stop[0]:
        break

  mapping_local: Dict[int, int] = {}
  total_bound = sum(peak_bound)
  backtrack(0, 0.0, total_bound, {})

  return ClusterResult(
    peaks=cluster_peaks,
    hypotheses=prune_hypotheses(hypotheses, total_peaks=total_peaks, score_factor=score_factor, max_keep=max_keep),
  )


def _compatible_with_usage(
  base_mapping: Dict[int, int],
  base_usage: Counter,
  extra_mapping: Dict[int, int],
  site_capacities: Dict[int, int],
) -> bool:
  """Compatibility check reusing a base's precomputed site-usage counter."""

  added = Counter()
  for peak_index, site_index in extra_mapping.items():
    if peak_index in base_mapping:
      if base_mapping[peak_index] != site_index:
        return False
      continue
    added[site_index] += 1
    if base_usage[site_index] + added[site_index] > site_capacities[site_index]:
      return False
  return True


def mappings_are_compatible(
  base_mapping: Dict[int, int],
  extra_mapping: Dict[int, int],
  site_capacities: Dict[int, int],
) -> bool:
  return _compatible_with_usage(base_mapping, Counter(base_mapping.values()), extra_mapping, site_capacities)


def _single_extra_peak(extra_hypotheses: Sequence[AssignmentHypothesis]):
  """Return the shared peak index if every extra maps exactly one, same peak."""

  peak = None
  for extra_hypothesis in extra_hypotheses:
    if len(extra_hypothesis.mapping) != 1:
      return None
    only_peak = next(iter(extra_hypothesis.mapping))
    if peak is None:
      peak = only_peak
    elif peak != only_peak:
      return None
  return peak


def _materialize_orphan(base_hypothesis, orphan_peak, site, score):
  mapping = dict(base_hypothesis.mapping)
  if site is not None:
    mapping[orphan_peak] = site
  return AssignmentHypothesis(mapping=mapping, score=score)


def _merge_single_peak(
  base_hypotheses: Sequence[AssignmentHypothesis],
  extra_hypotheses: Sequence[AssignmentHypothesis],
  orphan_peak: int,
  site_capacities: Dict[int, int],
  score_matrix: np.ndarray,
  model_matrix: np.ndarray,
  total_peaks: int,
  score_factor: float,
  max_keep: int,
) -> List[AssignmentHypothesis]:
  """Extend every base hypothesis by a single orphan peak.

  Candidate scores are computed first (the orphan only interacts with its few
  non-zero NOE partners, so each base needs a handful of terms, not all ~k),
  the keep-band is applied to the *scores*, and only the surviving candidates
  are materialised into full mappings.  Enumeration built ~25x more
  hypotheses than survive pruning; deferring the dict copies is the main win.
  Output order and values match the generic merge path.
  """

  score_row = score_matrix[orphan_peak]
  score_col = score_matrix[:, orphan_peak]
  self_term = float(score_matrix[orphan_peak, orphan_peak])

  partners = set(np.nonzero(score_row)[0].tolist()) | set(np.nonzero(score_col)[0].tolist())
  partners.discard(orphan_peak)
  partner_list = sorted(partners)
  partner_row = {peak: float(score_row[peak]) for peak in partner_list}
  partner_col = {peak: float(score_col[peak]) for peak in partner_list}

  unique_sites = sorted({extra_hypothesis.mapping[orphan_peak] for extra_hypothesis in extra_hypotheses})
  site_array = np.array(unique_sites, dtype=int)
  self_delta = self_term * model_matrix[site_array, site_array]

  # candidates: (score, mapping_length, base_hypothesis, site or None)
  candidates = []
  has_orphan_in_base = False
  for base_hypothesis in base_hypotheses:
    base_mapping = base_hypothesis.mapping
    if orphan_peak in base_mapping:
      # Never observed in practice (orphans are unassigned peaks); handle it
      # correctly via the safe dedup path below rather than the fast band cut.
      has_orphan_in_base = True
      fixed_site = base_mapping[orphan_peak]
      for extra_hypothesis in extra_hypotheses:
        if extra_hypothesis.mapping[orphan_peak] == fixed_site:
          candidates.append((base_hypothesis.score, len(base_mapping), base_hypothesis, None))
      continue

    rel_sites = []
    rel_row = []
    rel_col = []
    for peak in partner_list:
      site = base_mapping.get(peak)
      if site is not None:
        rel_sites.append(site)
        rel_row.append(partner_row[peak])
        rel_col.append(partner_col[peak])

    delta = self_delta
    if rel_sites:
      rel_array = np.array(rel_sites, dtype=int)
      delta = (
        delta
        + model_matrix[np.ix_(site_array, rel_array)] @ np.array(rel_row)
        + model_matrix[np.ix_(rel_array, site_array)].T @ np.array(rel_col)
      )
    score_by_site = {site: base_hypothesis.score + float(delta[i]) for i, site in enumerate(unique_sites)}

    length = len(base_mapping) + 1
    usage = Counter(base_mapping.values())
    for extra_hypothesis in extra_hypotheses:
      site = extra_hypothesis.mapping[orphan_peak]
      if usage[site] >= site_capacities[site]:
        continue
      candidates.append((score_by_site[site], length, base_hypothesis, site))

  if not candidates:
    return []
  if has_orphan_in_base:
    # Distinct mappings are not guaranteed; fall back to materialise-then-prune.
    hypotheses = [_materialize_orphan(base, orphan_peak, site, score) for score, _, base, site in candidates]
    return prune_hypotheses(hypotheses, total_peaks=total_peaks, score_factor=score_factor, max_keep=max_keep)

  # Fast path: all combined mappings are distinct, so pruning needs only the
  # scores.  Mirror prune_hypotheses' cutoff/sort exactly, then materialise.
  best_score = max(score for score, _, _, _ in candidates)
  assigned_ratio = max(length for _, length, _, _ in candidates) / max(1, total_peaks)
  slack = relative_score_slack(score_factor, assigned_ratio)
  cutoff = best_score * (1.0 - slack) if best_score > 0 else best_score - max(1.0, slack)
  survivors = [candidate for candidate in candidates if candidate[0] >= cutoff]
  survivors.sort(key=lambda candidate: (-candidate[0], candidate[1]))
  return [_materialize_orphan(base, orphan_peak, site, score) for score, _, base, site in survivors[:max_keep]]


def merge_hypothesis_sets(
  base_hypotheses: Sequence[AssignmentHypothesis],
  extra_hypotheses: Sequence[AssignmentHypothesis],
  site_capacities: Dict[int, int],
  score_matrix: np.ndarray,
  model_matrix: np.ndarray,
  total_peaks: int,
  score_factor: float,
  max_keep: int = 512,
) -> List[AssignmentHypothesis]:
  orphan_peak = _single_extra_peak(extra_hypotheses) if extra_hypotheses else None
  if orphan_peak is not None:
    return _merge_single_peak(
      base_hypotheses, extra_hypotheses, orphan_peak, site_capacities,
      score_matrix, model_matrix, total_peaks, score_factor, max_keep,
    )

  merged: List[AssignmentHypothesis] = []
  for base_hypothesis in base_hypotheses:
    base_mapping = base_hypothesis.mapping
    base_peaks = np.fromiter(base_mapping.keys(), dtype=int, count=len(base_mapping))
    base_sites = np.fromiter(base_mapping.values(), dtype=int, count=len(base_mapping))
    base_usage = Counter(base_mapping.values())
    for extra_hypothesis in extra_hypotheses:
      if not _compatible_with_usage(base_mapping, base_usage, extra_hypothesis.mapping, site_capacities):
        continue
      new_items = [(peak, site) for peak, site in extra_hypothesis.mapping.items() if peak not in base_mapping]
      combined_mapping = dict(base_mapping)
      combined_mapping.update(extra_hypothesis.mapping)
      # Combined score = base score + contributions of the newly-added peaks
      # (their self term, cross terms with the base, and among themselves),
      # avoiding an O(k^2) rescore of the whole mapping.
      delta = 0.0
      cur_peaks, cur_sites = base_peaks, base_sites
      for peak, site in new_items:
        delta += float(score_matrix[peak, peak] * model_matrix[site, site])
        if cur_peaks.size:
          delta += float(score_matrix[peak, cur_peaks] @ model_matrix[site, cur_sites])
          delta += float(score_matrix[cur_peaks, peak] @ model_matrix[cur_sites, site])
        cur_peaks = np.append(cur_peaks, peak)
        cur_sites = np.append(cur_sites, site)
      merged.append(AssignmentHypothesis(mapping=combined_mapping, score=base_hypothesis.score + delta))
  return prune_hypotheses(merged, total_peaks=total_peaks, score_factor=score_factor, max_keep=max_keep)


def connected_orphan_order(score_matrix: np.ndarray, assigned_peaks: Iterable[int], orphan_peaks: Iterable[int]) -> List[int]:
  assigned_list = list(assigned_peaks)
  orphan_list = list(orphan_peaks)
  if not assigned_list:
    return orphan_list
  ranking = []
  for peak_index in orphan_list:
    support = float(np.sum(score_matrix[peak_index, assigned_list]) + np.sum(score_matrix[assigned_list, peak_index]))
    ranking.append((support, peak_index))
  ranking.sort(reverse=True)
  return [peak_index for _, peak_index in ranking]
