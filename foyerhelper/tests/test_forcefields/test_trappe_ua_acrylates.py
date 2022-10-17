import unittest

from mbuildhelper.lib.recipes import PolyAlkylAcrylateUA
import foyerhelper


class TestTraPPEUAAcrylates(unittest.TestCase):

    def test_butyl_acrylate_3mer(self):
        pba_3mer = PolyAlkylAcrylateUA("3*A4")
        forcefield = foyerhelper.ForceField(forcefield_files="trappe-ua-acrylates.xml")
        pba_3mer = forcefield.apply(pba_3mer)

    def test_methyl_methacrylate_3mer(self):
        pba_3mer = PolyAlkylAcrylateUA("3*mA1")
        forcefield = foyerhelper.ForceField(forcefield_files="trappe-ua-acrylates.xml")
        pba_3mer = forcefield.apply(pba_3mer)

    def test_lauryl_methacrylate_5mer(self):
        plma_5mer = PolyAlkylAcrylateUA("5*mA12")
        forcefield = foyerhelper.ForceField(forcefield_files="trappe-ua-acrylates.xml")
        plma_5mer = forcefield.apply(plma_5mer)