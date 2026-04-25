"""Backward-compatible wrapper.

The canonical prep now lives at the project root in 06prep_disc.py
and writes size-tiered all_disc{,10,100,1000}.tsv files into
03data_bbknn_b/. This wrapper invokes it and copies the binary
subsets into the current directory so existing run.sh / run2.sh
keep working.
"""
import os
import runpy
import shutil
import sys
from pathlib import Path

if __name__ == "__main__":
    root = Path(__file__).resolve().parent.parent
    sys.argv = [sys.argv[0]]
    runpy.run_path(str(root / "06prep_disc.py"), run_name="__main__")

    src = root / "03data_bbknn_b"
    here = Path(__file__).resolve().parent
    for name in ("all_disc.tsv", "all_disc10.tsv", "all_disc100.tsv", "all_disc1000.tsv"):
        s = src / name
        if s.exists():
            shutil.copy2(s, here / name)
            print(f"copied -> {here / name}")
