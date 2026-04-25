"""Backward-compatible wrapper.

The canonical prep now lives at the project root in 06prep_disc.py
and writes size-tiered all_disc{,10,100,1000}.tsv (and the _tri
variants) into 03data_bbknn_b_<source>_<level>/, plus per-tissue
files under tissue[_tri]/. This wrapper invokes it for the canonical
discretization branch (--source r --level batch) and stages copies
into ./data_bin/ and ./data_tri/ so existing scripts keep working.
"""
import runpy
import shutil
import sys
from pathlib import Path

SOURCE = "r"
LEVEL  = "batch"


def _stage(src_root: Path, dest: Path, suffix: str, tissue_subdir: str):
    dest.mkdir(parents=True, exist_ok=True)
    for n in ("", "10", "100", "1000"):
        name = f"all_disc{suffix}{n}.tsv"
        s = src_root / name
        if s.exists():
            shutil.copy2(s, dest / name)
    src_tissue = src_root / tissue_subdir
    if src_tissue.is_dir():
        (dest / tissue_subdir).mkdir(exist_ok=True)
        for f in src_tissue.glob("*.tsv"):
            shutil.copy2(f, dest / tissue_subdir / f.name)


if __name__ == "__main__":
    root = Path(__file__).resolve().parent.parent
    sys.argv = [sys.argv[0], "--source", SOURCE, "--level", LEVEL]
    runpy.run_path(str(root / "06prep_disc.py"), run_name="__main__")

    src = root / f"03data_bbknn_b_{SOURCE}_{LEVEL}"
    here = Path(__file__).resolve().parent
    _stage(src, here / "data_bin", suffix="",     tissue_subdir="tissue")
    _stage(src, here / "data_tri", suffix="_tri", tissue_subdir="tissue_tri")
