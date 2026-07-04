"""Enable ``python3 -m magic ...`` as an entrypoint."""

from .cli import main


if __name__ == '__main__':
  raise SystemExit(main())
