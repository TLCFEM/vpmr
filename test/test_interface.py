import importlib
import sys
import unittest
from unittest import mock

import matplotlib
import numpy as np

matplotlib.use("Agg")


class PythonInterfaceTest(unittest.TestCase):
    def setUp(self):
        self.backend = mock.Mock()
        self.backend.vpmr.return_value = (
            [1 + 0j, 2 + 0j],
            [3 + 0j, 4 + 0j],
        )
        sys.modules["_pyvpmr"] = self.backend
        sys.modules.pop("pyvpmr", None)
        self.pyvpmr = importlib.import_module("pyvpmr")

    def tearDown(self):
        sys.modules.pop("pyvpmr", None)
        sys.modules.pop("_pyvpmr", None)

    def test_vpmr_accepts_named_kwargs(self):
        result = self.pyvpmr.vpmr(
            terms=12,
            max_exponent=6,
            precision_bits=320,
            quadrature_order=800,
            precision_multiplier=1.2,
            tolerance=1e-10,
            kernel="exp(-t)",
            omit_trivial_terms=False,
        )

        self.backend.vpmr.assert_called_once_with(
            n=12,
            c=6,
            d=320,
            q=800,
            m=1.2,
            e=1e-10,
            k="exp(-t)",
            omit=False,
        )
        self.assertIsInstance(result, self.pyvpmr.VPMRResult)
        self.assertTrue(np.array_equal(result.weights, np.array([1 + 0j, 2 + 0j])))
        self.assertTrue(np.array_equal(result.poles, np.array([3 + 0j, 4 + 0j])))

    def test_vpmr_accepts_options_object(self):
        options = self.pyvpmr.VPMROptions(terms=18, kernel="exp(-t)")
        self.pyvpmr.vpmr(options)
        self.backend.vpmr.assert_called_once_with(
            n=18,
            c=4,
            d=0,
            q=500,
            m=1.05,
            e=1e-8,
            k="exp(-t)",
            omit=True,
        )

    def test_vpmr_rejects_options_and_kwargs_together(self):
        with self.assertRaisesRegex(ValueError, "either 'options' or keyword"):
            self.pyvpmr.vpmr(self.pyvpmr.VPMROptions(), terms=10)

    def test_result_helpers_accept_result_object(self):
        result = self.pyvpmr.VPMRResult([2.0 + 0j], [0.5 + 0j])

        self.assertAlmostEqual(result.evaluate(2.0).real, 2.0 * np.exp(-1.0))
        values = result.evaluate(np.array([0.0, 1.0]))
        self.assertEqual(values.shape, (2,))
        self.assertIn("NonviscousNewmark", result.to_global_damping())
        self.assertIn("ElementalNonviscous", result.to_elemental_damping())

    def test_options_validation(self):
        with self.assertRaisesRegex(ValueError, "positive integer"):
            self.pyvpmr.VPMROptions(terms=0)
        with self.assertRaisesRegex(ValueError, "positive finite number"):
            self.pyvpmr.VPMROptions(tolerance=0.0)

    def test_plot_accepts_result_object(self):
        result = self.pyvpmr.VPMRResult([1.0 + 0j], [1.0 + 0j])

        fig, (ax1, ax2) = self.pyvpmr.plot(
            result,
            lambda x: np.exp(-x),
            show=False,
        )

        self.assertEqual(ax1.get_ylabel(), "kernel function $g(t)$")
        self.assertEqual(ax2.get_ylabel(), "absolute error")
        fig.clf()

    def test_helper_functions_require_result(self):
        with self.assertRaisesRegex(TypeError, "VPMRResult"):
            self.pyvpmr.to_global_damping(([1 + 0j], [2 + 0j]))


if __name__ == "__main__":
    unittest.main()
