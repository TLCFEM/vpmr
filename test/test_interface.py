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

    def test_descriptive_keywords_return_named_result(self):
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
        weights, poles = result
        self.assertTrue(np.array_equal(weights, result.m))
        self.assertTrue(np.array_equal(poles, result.s))

    def test_duplicate_keyword_aliases_must_match(self):
        with self.assertRaisesRegex(ValueError, "must match"):
            self.pyvpmr.vpmr(n=3, terms=4)

    def test_result_helpers_accept_result_object(self):
        result = self.pyvpmr.VPMRResult([2.0 + 0j], [0.5 + 0j])

        self.assertAlmostEqual(result.evaluate(2.0).real, 2.0 * np.exp(-1.0))
        values = result.evaluate(np.array([0.0, 1.0]))
        self.assertEqual(values.shape, (2,))
        self.assertIn("NonviscousNewmark", result.to_global_damping())
        self.assertIn("ElementalNonviscous", result.to_elemental_damping())

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


if __name__ == "__main__":
    unittest.main()
