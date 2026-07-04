"""Command-line entrypoint for the clean MAGIC implementation.

Supports three subcommands (``run`` / ``validate`` / ``template``) while
preserving the original ``Magic3v1.py input [output_dir]`` calling style.
"""

import argparse
from pathlib import Path
import sys

from .config import emit_template
from .inputs import validate_magic_inputs
from .pipeline import run_magic


def build_cli_parser() -> argparse.ArgumentParser:
  parser = argparse.ArgumentParser(
    description='Run the MAGIC methyl-assignment workflow or inspect its input files.',
    epilog=(
      'Examples:\n'
      '  python3 Magic3v1.py run magic_control.txt --output-dir run_01\n'
      '  python3 Magic3v1.py validate magic_control.txt\n'
      '  python3 Magic3v1.py template magic_control.txt'
    ),
    formatter_class=argparse.RawDescriptionHelpFormatter,
  )
  subparsers = parser.add_subparsers(dest='command')

  run_parser = subparsers.add_parser('run', help='run the MAGIC pipeline from a control file')
  run_parser.add_argument('input_file', help='control file describing the HMQC/NOESY/PDB inputs')
  run_parser.add_argument(
    '--output-dir',
    help='output directory (relative paths are safest for this legacy codebase)',
  )

  validate_parser = subparsers.add_parser('validate', help='validate a control file and the files it references')
  validate_parser.add_argument('input_file', help='control file to validate')

  template_parser = subparsers.add_parser('template', help='print or write a starter control file')
  template_parser.add_argument(
    'output_file',
    nargs='?',
    help='optional path for writing the template instead of printing it',
  )
  return parser


def normalize_legacy_cli_args(argv):
  """Preserve the original `Magic3v1.py input [output_dir]` calling style."""

  if not argv:
    return argv
  if argv[0] in {'run', 'validate', 'template', '-h', '--help'}:
    return argv
  normalized = ['run', argv[0]]
  if len(argv) > 1:
    normalized.extend(['--output-dir', argv[1]])
  if len(argv) > 2:
    normalized.extend(argv[2:])
  return normalized


def parse_cli_args(argv):
  parser = build_cli_parser()
  args = parser.parse_args(normalize_legacy_cli_args(argv))
  if args.command is None:
    parser.print_help()
    parser.exit(2)
  return args


def resolve_control_file(input_file) -> Path:
  control_file = Path(input_file).expanduser()
  control_file = control_file.resolve() if control_file.is_absolute() else (Path.cwd() / control_file).resolve()
  if not control_file.exists():
    raise SystemExit(f'Error: Control file not found: {control_file}')
  return control_file


def main(argv=None):
  args = parse_cli_args(sys.argv[1:] if argv is None else argv)
  if args.command == 'template':
    emit_template(args.output_file)
    return 0

  control_file = resolve_control_file(args.input_file)
  if args.command == 'validate':
    validate_magic_inputs(control_file)
    return 0

  run_magic(control_file, args.output_dir)
  return 0
