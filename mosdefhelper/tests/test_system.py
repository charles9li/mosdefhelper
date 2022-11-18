import unittest

import mbuild as mb

from mbuildhelper.lib.molecules import DodecaneUA
from mbuildhelper.lib.recipes import PolyAlkylAcrylateUA
import foyerhelper
import mosdefhelper


class TestSystem(unittest.TestCase):

    def test_pba_5mer_dodecane(self):
        # create system
        box = mb.Box([5, 5, 5])
        system = mosdefhelper.System(box=box)

        # add compounds and corresponding force fields
        print("adding compounds")
        forcefield = foyerhelper.Forcefield(name='trappe-ua-acrylates')
        system.add_compound(PolyAlkylAcrylateUA("5*A4"), n=2,
                            forcefield=forcefield, apply_kwargs={'residues': "A4"})
        system.add_compound(DodecaneUA(), n=2,
                            forcefield=forcefield, apply_kwargs={'residues': "DodecaneUA"})

        # convert to parmed and save
        print("converting to parmed")
        system_pmd = system.to_parmed()
        system_pmd.save("pba5mer_dod.gro", overwrite=True)
        system_pmd.save("pba5mer_dod.top", overwrite=True)
        system_pmd.save("pba5mer_dod.pdb", overwrite=True)


if __name__ == '__main__':
    unittest.main()
