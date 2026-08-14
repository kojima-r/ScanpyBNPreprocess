import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[1]


class PipelineIntegrationTest(unittest.TestCase):
    def test_pseudo_bulk_continuous_and_discrete_branches(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            inputs = base / "input"
            inputs.mkdir()
            (inputs / "Heart.txt").write_text(
                "@name\tG1\tG2\tG3\n"
                "Heart|3m|b1|c1\t0\t1\t5\n"
                "Heart|3m|b1|c2\t2\t3\t5\n"
                "Heart|18m|b2|c3\t4\t5\t5\n",
                encoding="utf-8",
            )
            (inputs / "Lung.txt").write_text(
                "@name\tG1\tG2\tG3\n"
                "Lung|3m|b1|c1\t1\t2\t5\n"
                "Lung|18m|b2|c2\t5\t6\t5\n",
                encoding="utf-8",
            )
            config = {
                "dataset": {"mode": "bbknn", "target": "toy"},
                "pipeline": {"preprocess": False, "aggregation": True, "continuous": True, "discrete": True},
                "aggregation": {
                    "method": "pseudo_bulk", "level": "age", "workers": 1,
                    "input_glob": "input/*.txt", "output_dir": "aggregated",
                },
                "continuous": {"transpose_output_dir": "transposed", "merge_output": "continuous/all.txt"},
                "discrete": {
                    "merge_output": "discrete/all.txt", "discretize_output_dir": "discrete",
                    "prepare_output_dir": "prepared",
                },
            }
            config_path = base / "config.yaml"
            config_path.write_text(yaml.safe_dump(config), encoding="utf-8")

            result = subprocess.run(
                [sys.executable, str(ROOT / "pipeline.py"), "--config", str(config_path)],
                cwd=ROOT, text=True, capture_output=True,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertTrue((base / "continuous/all.txt").exists())
            self.assertTrue((base / "prepared/all_disc.tsv").exists())
            binary = pd.read_csv(base / "prepared/all_disc.tsv", sep="\t")
            self.assertEqual(list(binary.columns), ["G1", "G2"])
            self.assertEqual(binary.shape, (4, 2))


if __name__ == "__main__":
    unittest.main()
