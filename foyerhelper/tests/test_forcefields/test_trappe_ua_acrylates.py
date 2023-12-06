import unittest

from mbuildhelper.lib.recipes import PolyAlkylAcrylateUA
import foyerhelper


class TestTraPPEUAAcrylates(unittest.TestCase):

    def test_butyl_acrylate_monomer(self):
        ba = PolyAlkylAcrylateUA("1*A4")
        forcefield = foyerhelper.Forcefield(name="trappe-ua-acrylates")
        ba = forcefield.apply(ba, infer_residues=True, assert_angle_params=True, assert_dihedral_params=True)
        self.assertListEqual([0.0, 0.0, 0.4, -0.4, -0.25, 0.25, 0.0, 0.0, 0.0], [a.charge for a in ba.atoms])

    def test_butyl_acrylate_3mer(self):
        pba_3mer = PolyAlkylAcrylateUA("3*A4")
        forcefield = foyerhelper.Forcefield(name="trappe-ua-acrylates")
        pba_3mer = forcefield.apply(pba_3mer, infer_residues=True)
        self.assertListEqual([0.0, 0.0, 0.4, -0.4, -0.25, 0.25, 0.0, 0.0, 0.0]*3, [a.charge for a in pba_3mer.atoms])

    def test_butyl_acrylate_5mer(self):
        pba_5mer = PolyAlkylAcrylateUA("5*A4")
        forcefield = foyerhelper.Forcefield(name="trappe-ua-acrylates")
        pba_5mer = forcefield.apply(pba_5mer, infer_residues=True)
        self.assertListEqual([0.0, 0.0, 0.4, -0.4, -0.25, 0.25, 0.0, 0.0, 0.0]*5, [a.charge for a in pba_5mer.atoms])

    def test_dodecyl_acrylate_monomer(self):
        dda = PolyAlkylAcrylateUA("1*A12")
        forcefield = foyerhelper.Forcefield(name="trappe-ua-acrylates")
        dda = forcefield.apply(dda, infer_residues=True, assert_angle_params=True, assert_dihedral_params=True)
        self.assertListEqual([0.0, 0.0, 0.4, -0.4, -0.25, 0.25] + [0.0]*11, [a.charge for a in dda.atoms])

    def test_methyl_methacrylate_3mer(self):
        pmma_3mer = PolyAlkylAcrylateUA("3*mA1")
        forcefield = foyerhelper.Forcefield(name="trappe-ua-acrylates")
        pmma_3mer = forcefield.apply(pmma_3mer, infer_residues=True)
        self.assertListEqual(["mA1"]*3, [r.name for r in pmma_3mer.residues])

    def test_lauryl_methacrylate_5mer(self):
        plma_5mer = PolyAlkylAcrylateUA("5*mA12")
        forcefield = foyerhelper.Forcefield(name="trappe-ua-acrylates")
        plma_5mer = forcefield.apply(plma_5mer, infer_residues=True)
        self.assertListEqual(["mA12"]*5, [r.name for r in plma_5mer.residues])


if __name__ == '__main__':
    unittest.main()
