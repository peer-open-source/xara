#!/usr/bin/env python
"""
./win_repair.py ./win32/ wheelhouse/*.whl

Claudio Perez
"""
import sys
import zipfile
from pathlib import Path

bindir = Path("opensees/bin/")

for build in Path(".").glob(sys.argv[2]):
    print(f"Repairing {build}", file=sys.stderr)
    with zipfile.ZipFile(build, 'a') as zipf:
        for file in Path(sys.argv[1]).glob("*.dll"):
            print(f" :: adding {file.name} to {bindir/file.name}", file=sys.stderr)
            zipf.write(file, bindir/file.name)

