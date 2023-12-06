import unittest

from mbuildhelper.lib.molecules import StearicAcidUA
import foyerhelper


class TestTraPPEUACarboxylicAcids(unittest.TestCase):

    def test_stearic_acid(self):
        stearic_acid = StearicAcidUA()
        forcefield = foyerhelper.Forcefield(name='trappe-ua-carboxylic-acids')
        stearic_acid = forcefield.apply(stearic_acid)
        self.assertListEqual([-0.46, 0.37, 0.42, -0.45, 0.12]+[0.0]*16, [a.charge for a in stearic_acid.atoms])


if __name__ == '__main__':
    unittest.main()
