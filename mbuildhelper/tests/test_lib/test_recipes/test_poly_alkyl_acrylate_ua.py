import unittest

from mbuildhelper.lib.recipes import PolyAlkylAcrylateUA


class TestPolyAlkylAcrylateUA(unittest.TestCase):

    def test_butyl_acrylate_3mer(self):
        pba_3mer_ua = PolyAlkylAcrylateUA("3*A4")
        pba_3mer_ua.save("pba_3mer_ua.mol2", overwrite=True)

    def test_butyl_acrylate_5mer(self):
        pba_5mer_ua = PolyAlkylAcrylateUA("5*A4")
        pba_5mer_ua.save("pba_5mer_ua.mol2", overwrite=True)

    def test_lauryl_acrylate_4mer(self):
        pla_4mer_ua = PolyAlkylAcrylateUA("4*A12")
        pla_4mer_ua.save("pla_4mer_ua.mol2", overwrite=True)

    def test_lauryl_methacrylate_4mer(self):
        plma_4mer_ua = PolyAlkylAcrylateUA("4*mA12")
        plma_4mer_ua.save("plma_4mer_ua.mol2", overwrite=True)


if __name__ == '__main__':
    unittest.main()
