import unittest

import warnings

import mbuild as mb

from mbuildhelper.lib.molecules import DodecaneUA
import foyerhelper


class TestDodecaneUA(unittest.TestCase):

    def test_dodecane_ua(self):
        d = DodecaneUA()
        self.assertEqual(12, d.n_particles)
        self.assertEqual(11, d.n_bonds)
        d.save("dodecane_ua.mol2", overwrite=True)

    def test_fill(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            box = mb.Box(lengths=[5, 5, 5])
            dodecane = mb.fill_box(compound=DodecaneUA(), density=750, box=box)
            trappeua = foyerhelper.Forcefield(name="trappe-ua-acrylates")
            dodecane = trappeua.apply(dodecane)
            dodecane.save("dodecane.gro", overwrite=True)
            dodecane.save("dodecane.top", overwrite=True)


if __name__ == '__main__':
    unittest.main()
