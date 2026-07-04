"""Resolving, validating and copying the files a control file points at."""

import os
from pathlib import Path

from .config import MagicControlFile, parse_magic_control_file


def resolve_input_path(input_name: str, base_directory: Path) -> Path:
  path = Path(input_name).expanduser()
  if not path.is_absolute():
    path = base_directory / path
  return path


def resolve_magic_inputs(control_file: Path, config: MagicControlFile):
  base_directory = control_file.parent
  resolved_inputs = {
    'control': control_file,
    'hmqc': resolve_input_path(config.hmqc_name, base_directory),
    'pdb': resolve_input_path(config.pdb_name, base_directory),
    'seq': resolve_input_path(config.seq_name, base_directory),
    'noesy': [resolve_input_path(item, base_directory) for item in config.noesy_list],
  }

  missing_paths = []
  for key, value in resolved_inputs.items():
    if key == 'noesy':
      missing_paths.extend(str(path) for path in value if not path.exists())
    elif not value.exists():
      missing_paths.append(str(value))
  if missing_paths:
    raise FileNotFoundError(
      'Missing input files referenced by {}:\n- {}'.format(
        control_file,
        '\n- '.join(missing_paths),
      )
    )
  return resolved_inputs


def resolve_control_inputs(control_file: Path):
  config = parse_magic_control_file(control_file)
  resolved = resolve_magic_inputs(control_file, config)
  return config, resolved


def validate_magic_inputs(control_file: Path):
  config = parse_magic_control_file(control_file)
  resolved_inputs = resolve_magic_inputs(control_file, config)
  print('Control file OK: {}'.format(control_file))
  print('HMQC: {}'.format(resolved_inputs['hmqc']))
  print('PDB: {}'.format(resolved_inputs['pdb']))
  print('SEQ: {}'.format(resolved_inputs['seq']))
  for noesy_path in resolved_inputs['noesy']:
    print('NOESY: {}'.format(noesy_path))


def default_output_directory(start_time) -> str:
  return 'magic_{}'.format(start_time.strftime('%Y%m%d_%H%M%S'))


def normalize_output_directory_name(output_dir, start_time) -> str:
  raw_output = output_dir if output_dir else default_output_directory(start_time)
  output_path = Path(raw_output).expanduser()
  if output_path.is_absolute():
    output_path = Path(os.path.relpath(output_path, Path.cwd()))
  sanitized_parts = [
    part if part in {'.', '..'} else part.replace('.', '_')
    for part in output_path.parts
  ]
  return os.path.join(*sanitized_parts)
