import unittest

import mbuild as mb
import gmso
from gmso.external.convert_parmed import from_parmed
import parmed as pmd

from mbuildhelper.lib.surfaces import Hematite0001Surface
import foyerhelper


class TestIronOxide(unittest.TestCase):

    def test_hematite(self):
        hematite: mb.Compound = Hematite0001Surface(lengths=[2, 2, 1])
        forcefield = foyerhelper.Forcefield(name="iron-oxide")
        hematite: pmd.Structure = forcefield.apply(hematite,
                                                   assert_angle_params=False,
                                                   assert_dihedral_params=False,
                                                   set_masses_to_zero=True)
        hematite.save("hematite.top", overwrite=True)
        hematite.save("hematite.gro", overwrite=True)
        hematite: gmso.Topology = from_parmed(hematite)
        for bt in hematite.bond_types:
            print(bt.member_types)


if __name__ == '__main__':
    unittest.main()
