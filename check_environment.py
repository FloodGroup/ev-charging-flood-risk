"""Check the Python packages used to draw figures and run the analyses."""

import argparse
from importlib import import_module
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
import sys

from shared import paths


REQUIREMENTS = Path(__file__).with_name("requirements.txt")
IMPORT_NAMES = {"scikit-learn": "sklearn", "pillow": "PIL"}


def check_packages():
    """Report installed and pinned versions. Return False for import or lookup failures."""
    available = True
    for line in REQUIREMENTS.read_text(encoding="utf-8").splitlines():
        requirement = line.partition("#")[0].strip()
        if not requirement:
            continue
        name, pinned = requirement.split("==")
        try:
            import_module(IMPORT_NAMES.get(name, name))
            installed = version(name)
        except (ImportError, OSError, PackageNotFoundError) as error:
            print(f"[missing] {name}  {error}")
            available = False
            continue
        status = "ok" if installed == pinned else f"expected {pinned}"
        print(f"[{status}] {name} {installed}")
    return available


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", action="store_true", help="list local analysis inputs and outputs")
    args = parser.parse_args()

    print(f"Python {sys.version.split()[0]}  ({sys.executable})")
    supported = sys.version_info >= (3, 11)
    if not supported:
        print("Python 3.11 or later is required.")
    available = check_packages()
    if args.inputs:
        for key, path, present in paths.describe():
            status = "present" if present else "absent"
            print(f"[{status}] {key}  {path}")
    return 0 if supported and available else 1


if __name__ == "__main__":
    raise SystemExit(main())
