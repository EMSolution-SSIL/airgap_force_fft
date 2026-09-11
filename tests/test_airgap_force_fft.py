from __future__ import annotations

import tempfile
import unittest
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from airgap_fft import (
    dominant_modes,
    extract_torque_references,
    fft2_spectrum,
    integrated_quantities,
)
from airgap_io import (
    expand_sector,
    load_airgap_csv,
    load_airgap_json,
)


class AirgapForceFftTests(unittest.TestCase):
    def test_sector_expansion_makes_full_turn(self) -> None:
        theta = np.array([15.0, 45.0])
        br = np.array([[1.0, 2.0]])
        bt = np.array([[0.1, 0.2]])
        theta_full, br_full, bt_full = expand_sector(theta, br, bt, 60.0, 6)

        self.assertEqual(len(theta_full), 12)
        self.assertAlmostEqual(theta_full[0], 15.0)
        self.assertAlmostEqual(theta_full[-1], 345.0)
        np.testing.assert_allclose(br_full[0, :4], [1.0, 2.0, 1.0, 2.0])
        np.testing.assert_allclose(bt_full[0, :4], [0.1, 0.2, 0.1, 0.2])

    def test_anti_periodic_sector_expansion_alternates_sign(self) -> None:
        theta = np.array([15.0, 45.0])
        br = np.array([[1.0, 2.0]])
        bt = np.array([[0.1, 0.2]])
        theta_full, br_full, bt_full = expand_sector(
            theta, br, bt, 60.0, 6, symmetricity=1
        )

        self.assertEqual(len(theta_full), 12)
        np.testing.assert_allclose(br_full[0, :6], [1.0, 2.0, -1.0, -2.0, 1.0, 2.0])
        np.testing.assert_allclose(bt_full[0, :6], [0.1, 0.2, -0.1, -0.2, 0.1, 0.2])

    def test_uniform_tangential_stress_torque(self) -> None:
        theta = np.linspace(0.0, 2.0 * np.pi, 8, endpoint=False)
        sigma_r = np.zeros((3, len(theta)))
        sigma_t = np.full((3, len(theta)), 10.0)
        metrics = integrated_quantities(theta, sigma_r, sigma_t, radius=2.0, stack_length=3.0)

        np.testing.assert_allclose(metrics["Torque_Nm"], 10.0 * 2.0**2 * 3.0 * 2.0 * np.pi)
        np.testing.assert_allclose(metrics["UMP_N"], 0.0)

    def test_traveling_wave_peak_uses_signed_axes(self) -> None:
        nt = 32
        nth = 48
        step = np.arange(nt, dtype=float)
        theta = np.linspace(0.0, 2.0 * np.pi, nth, endpoint=False)
        field = np.cos(3.0 * theta[None, :] - 2.0 * np.pi * 5.0 * step[:, None] / nt)

        sf, sm, coeff = fft2_spectrum(field, step, theta)
        modes = dominant_modes(sf, sm, coeff, n=2)
        pairs = set(zip(np.round(modes["step_frequency_Hz_or_order"], 8), modes["spatial_order_m"]))

        self.assertIn((-5.0 / nt, 3.0), pairs)
        self.assertIn((5.0 / nt, -3.0), pairs)

    def test_load_emsolution_json_gl80_shape(self) -> None:
        path = (
            Path(__file__).resolve().parents[1]
            / "data"
            / "GL80"
            / "healthy_position"
            / "motorGapB.json"
        )
        if not path.exists():
            self.skipTest("GL80 sample JSON not present")

        with path.open("r", encoding="utf-8") as f:
            raw = json.load(f)
        data = load_airgap_json(path, sector_periods="auto", resample_spatial=True)

        self.assertEqual(data.br.shape[0], raw["timeStep"]["numSteps"])
        self.assertEqual(data.bt.shape, data.br.shape)
        self.assertEqual(data.metadata["sector_periods"], 1)
        original_theta_count = len(raw["postData"]["gapB"]["position"]["theta"])
        self.assertEqual(data.br.shape[1], original_theta_count)
        expected_radius = np.mean(raw["postData"]["gapB"]["position"]["r"])
        self.assertAlmostEqual(data.radius or 0.0, expected_radius)
        self.assertAlmostEqual(data.stack_length or 0.0, 0.01)

    def test_load_csv_sector_shape(self) -> None:
        rows = [
            "time,theta_deg,Br,Btheta",
            "0.0,0.0,1.0,0.1",
            "0.0,30.0,2.0,0.2",
            "0.1,0.0,1.5,0.15",
            "0.1,30.0,2.5,0.25",
        ]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "sector.csv"
            path.write_text("\n".join(rows) + "\n", encoding="utf-8")
            data = load_airgap_csv(
                path,
                sector_periods=6,
                sector_angle_deg=60.0,
                symmetricity=1,
                resample_spatial=False,
            )

        self.assertEqual(data.br.shape, (2, 12))
        self.assertEqual(data.bt.shape, data.br.shape)
        self.assertEqual(data.metadata["sector_periods"], 6)
        self.assertEqual(data.metadata["symmetricity"], 1)

    def test_extract_pyemsol_force_nodal_torque_references(self) -> None:
        payload = {
            "timeStep": {"time": [0.0, 0.1, 0.2]},
            "postData": {
                "forceNodal": {
                    "forceNodalData": [
                        {"propertyNum": 10, "forceMZ": [1.0, 2.0, 3.0]},
                        {"propertyNum": "stator", "forceMZ": [-1.0, -2.0, -3.0]},
                        {"propertyNum": "rotor", "forceMZ": [1.1, 2.1, 3.1]},
                    ]
                }
            },
        }
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "output_transient.json"
            path.write_text(json.dumps(payload), encoding="utf-8")

            refs = extract_torque_references(path)

        self.assertEqual(list(refs.columns), [
            "time",
            "reference_Torque_Stator_Nm",
            "reference_Torque_Rotor_Nm",
        ])
        np.testing.assert_allclose(refs["time"], [0.0, 0.1, 0.2])
        np.testing.assert_allclose(refs["reference_Torque_Stator_Nm"], [-1.0, -2.0, -3.0])
        np.testing.assert_allclose(refs["reference_Torque_Rotor_Nm"], [1.1, 2.1, 3.1])


if __name__ == "__main__":
    unittest.main()
