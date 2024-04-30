import unittest

from mbuildhelper.lib.molecules import MethanolUA
import foyerhelper


class TestTraPPEUAAlcohols(unittest.TestCase):

    def test_methanol(self):
        methanol = MethanolUA()
        forcefield = foyerhelper.Forcefield(name='trappe-ua-alcohols')
        methanol = forcefield.apply(methanol)
        self.assertListEqual([-0.7, 0.435, 0.265], [a.charge for a in methanol.atoms])


if __name__ == '__main__':
    unittest.main()
