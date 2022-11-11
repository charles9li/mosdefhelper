import unittest

import mbuild as mb

from mbuildhelper.lib.recipes import PolyAlkylAcrylateUA
import foyerhelper


class TestPolyAlkylAcrylateUA(unittest.TestCase):

    def test_butyl_acrylate_3mer(self):
        pba_3mer_ua = PolyAlkylAcrylateUA("3*A4")
        pba_3mer_ua.save("pba_3mer_ua.mol2", overwrite=True)

    def test_butyl_acrylate_5mer(self):
        pba_5mer_ua = PolyAlkylAcrylateUA("5*A4")
        pba_5mer_ua.save("pba_5mer_ua.mol2", overwrite=True)

    def test_butyl_acrylate_10mer(self):
        pba_10mer_ua = PolyAlkylAcrylateUA("10*A4")
        box = mb.Box(lengths=[4, 4, 4])
        system = mb.fill_box(pba_10mer_ua, n_compounds=1, box=box)
        trappeua = foyerhelper.ForceField(forcefield_files="trappe-ua-acrylates.xml")
        system = trappeua.apply(system)
        system.save("pba_10mer.gro", overwrite=True)
        system.save("pba_10mer.top", overwrite=True)

    def test_lauryl_acrylate_4mer(self):
        pla_4mer_ua = PolyAlkylAcrylateUA("4*A12")
        pla_4mer_ua.save("pla_4mer_ua.mol2", overwrite=True)

    def test_lauryl_methacrylate_4mer(self):
        plma_4mer_ua = PolyAlkylAcrylateUA("4*mA12")
        plma_4mer_ua.save("plma_4mer_ua.mol2", overwrite=True)


if __name__ == '__main__':
    unittest.main()
