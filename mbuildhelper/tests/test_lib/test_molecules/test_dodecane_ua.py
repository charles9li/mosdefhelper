import unittest

from mbuildhelper.lib.molecules import DodecaneUA


class TestDodecaneUA(unittest.TestCase):

    def test_dodecane_ua(self):
        d = DodecaneUA()
        self.assertEqual(12, d.n_particles)
        self.assertEqual(11, d.n_bonds)
        d.save("dodecane_ua.mol2", overwrite=True)


if __name__ == '__main__':
    unittest.main()
