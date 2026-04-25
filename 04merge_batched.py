"""Backward-compatible wrapper: merge with --batched."""
import runpy
import sys

if __name__ == "__main__":
    sys.argv = [sys.argv[0], "--batched"] + sys.argv[1:]
    runpy.run_path("04merge.py", run_name="__main__")
