"""Clean-room MAGIC implementation based on Monneau et al. (2017)."""

from .cli import main
from .pipeline import run_magic

__all__ = ['main', 'run_magic']
