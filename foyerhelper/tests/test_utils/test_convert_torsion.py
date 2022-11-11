import unittest

import numpy as np

from foyerhelper.utils import convert_torsion


class TestConvertTorsion(unittest.TestCase):

    def test_trappeua(self):
        # opls coefficients for CHx-[CH2]-[CH2]-CHy torsion in kB units
        coefficients = np.array([0.0, 335.03, -68.19, 791.32])

        # convert to power form in kcal/mol units
        coefficients_power_kcalpermol = convert_torsion(coefficients, 'opls', 'power', 'kB', 'kcal/mol')
        self.assertIsNone(np.testing.assert_array_almost_equal(coefficients_power_kcalpermol,
                                                               [1.967, -4.052, 0.271, 6.290],
                                                               decimal=3))

        # convert to rb form in kJ/mol units
        coefficients_rb_kJpermol = convert_torsion(coefficients, 'opls', 'rb', 'kB', 'kJ/mol')
        self.assertIsNone(np.testing.assert_array_almost_equal(coefficients_rb_kJpermol,
                                                               [8.231068558, 16.95260727, 1.133926412, -26.31760224],
                                                               decimal=6))

        # now convert from power and kcal/mol units to rb and kJ/mol units to check consistency
        coefficients_rb_kJpermol_2 = convert_torsion(coefficients_power_kcalpermol, 'power', 'rb', 'kcal/mol', 'kJ/mol')
        self.assertIsNone(np.testing.assert_array_almost_equal(coefficients_rb_kJpermol_2, coefficients_rb_kJpermol))


if __name__ == '__main__':
    unittest.main()
