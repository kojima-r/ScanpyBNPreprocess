import tempfile
import time
import unittest
from pathlib import Path

import yaml

from server import RunManager


class RunManagerIntegrationTest(unittest.TestCase):
    def test_run_history_progress_log_and_results(self):
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
            config = {
                "dataset": {"family": "tms", "mode": "bbknn", "target": "ui_test"},
                "pipeline": {"preprocess": False, "aggregation": True, "discrete": True},
                "aggregation": {
                    "method": "pseudo_bulk",
                    "level": "age",
                    "workers": 1,
                    "input_glob": "input/*.txt",
                    "output_dir": "aggregated",
                },
                "discrete": {
                    "merge_output": "discrete/all.txt",
                    "discretize_output_dir": "discrete",
                    "prepare_output_dir": "prepared",
                },
            }
            config_path = base / "config.yaml"
            config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
            manager = RunManager(config_path, base / "runs")

            started = manager.start()
            deadline = time.monotonic() + 30
            while manager.status()["running"] and time.monotonic() < deadline:
                time.sleep(0.1)

            detail = manager.detail(started["id"])
            self.assertEqual(detail["status"], "succeeded", detail.get("log"))
            self.assertEqual(detail["progress"], 100)
            self.assertTrue(all(stage["status"] == "succeeded" for stage in detail["stages"]))
            self.assertIn("Pipeline completed.", detail["log"])
            self.assertGreater(detail["result_count"], 0)
            self.assertTrue(detail["files"])
            self.assertTrue(manager.result_file(detail["id"], detail["files"][0]["id"]).is_file())
            self.assertEqual(manager.list_runs()[0]["id"], detail["id"])

            manager.delete(detail["id"])
            self.assertEqual(manager.list_runs(), [])


if __name__ == "__main__":
    unittest.main()
