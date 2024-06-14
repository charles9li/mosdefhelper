import unittest

from mbuildhelper.lib.molecules import TolueneUA7Site, TolueneUA10Site
import foyerhelper


class TestTraPPEUAAlkylbezene(unittest.TestCase):

    def test_toluene_7site(self):
        toluene_ua_7site = TolueneUA7Site()
        toluene_ua_7site.name = "TOL"
        forcefield = foyerhelper.Forcefield(name="trappe-ua-alkylbenzene")
        toluene_ua_7site = forcefield.apply(
            toluene_ua_7site,
            residues="TOL"
        )

    def test_toluene_10site(self):
        toluene_ua_10site = TolueneUA10Site()
        toluene_ua_10site.name = "TOL"
        forcefield = foyerhelper.Forcefield(name="trappe-ua-alkylbenzene-modified")
        toluene_ua_10site = forcefield.apply(
            toluene_ua_10site,
            assert_angle_params=False,
            assert_dihedral_params=False,
            residues="TOL"
        )
        self.assertListEqual(
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.42, -1.21, -1.21, 0.0],
            [a.charge for a in toluene_ua_10site.atoms]
        )


if __name__ == '__main__':
    unittest.main()
