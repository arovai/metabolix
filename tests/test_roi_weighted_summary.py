import csv
import json
import tempfile
import unittest
from pathlib import Path

from metabolix.discovery import MRSInput
from metabolix.output_layout import OutputLayout, bids_filename
from metabolix.pipeline import _write_roi_weighted_estimates


class RoiWeightedSummaryTests(unittest.TestCase):
    def test_includes_overlap_weighted_component_qc_metrics(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            item = MRSInput(
                path=Path("sub-01_mrs.nii.gz"),
                sidecar=None,
                entities={"sub": "01"},
                standard_bids=True,
            )
            layout = OutputLayout(Path(temp_dir), item)
            layout.roi_dir.mkdir(parents=True)
            qc_path = layout.branch_qc_table("original")
            qc_path.parent.mkdir(parents=True)

            overlap_rows = [
                {"i": 0, "j": 0, "k": 0, "target_overlap_pct": 50},
                {"i": 1, "j": 0, "k": 0, "target_overlap_pct": 25},
                {"i": 2, "j": 0, "k": 0, "target_overlap_pct": 10},
                {"i": 3, "j": 0, "k": 0, "target_overlap_pct": 80},
            ]
            overlap_path = layout.roi_dir / bids_filename(
                item, "metabolix-roi-overlap", "voxelqc", ".tsv"
            )
            with overlap_path.open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=overlap_rows[0], delimiter="\t")
                writer.writeheader()
                writer.writerows(overlap_rows)

            metrics = {
                "raw_fit_NAA": (10, 30, 50, 70),
                "internal_reference_ratio_NAA": (1, 3, 5, 7),
                "crlb_percent_NAA": (10, 30, 50, 70),
                "snr_NAA": (5, 9, 13, 17),
                "linewidth_fwhm_hz_NAA": (4, 8, 12, 16),
            }
            voxel_rows = []
            for index in range(4):
                row = {
                    "i": index,
                    "j": 0,
                    "k": 0,
                    "status": "excluded" if index == 3 else "fit_completed",
                }
                row.update({name: values[index] for name, values in metrics.items()})
                voxel_rows.append(row)
            with qc_path.open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=voxel_rows[0], delimiter="\t")
                writer.writeheader()
                writer.writerows(voxel_rows)

            output_path = _write_roi_weighted_estimates(
                layout, ["original"], {"roi": {"minimum_overlap_percent": 20}}
            )
            with output_path.open(newline="") as stream:
                rows = list(csv.DictReader(stream, delimiter="\t"))
            by_metric = {row["metric"]: row for row in rows}

            self.assertEqual(set(by_metric), set(metrics))
            self.assertAlmostEqual(
                float(by_metric["crlb_percent_NAA"]["overlap_weighted_voxel_estimate"]),
                50 / 3,
            )
            self.assertAlmostEqual(
                float(by_metric["snr_NAA"]["overlap_weighted_voxel_estimate"]),
                19 / 3,
            )
            self.assertAlmostEqual(
                float(by_metric["linewidth_fwhm_hz_NAA"]["overlap_weighted_voxel_estimate"]),
                16 / 3,
            )
            self.assertTrue(all(row["weighted_voxel_count"] == "2" for row in rows))

            metadata = json.loads(output_path.with_suffix(".json").read_text())
            self.assertIn("snr_NAA", metadata["Metrics"])


if __name__ == "__main__":
    unittest.main()
