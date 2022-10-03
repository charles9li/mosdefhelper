import unittest

import mbuild as mb

from mbuildhelper.lib.recipes import AlkaneUA


class TestAlkaneUA(unittest.TestCase):

    def test_methane(self):
        methane = AlkaneUA(n=1)
        self.assertEqual(1, methane.n_particles)
        self.assertEqual("_CH4", methane[0].name)
        self.assertEqual(0, len(methane.all_ports()))

    def test_n1_cap_front(self):
        ch3_cap_front = AlkaneUA(n=1, cap_front=True, cap_end=False)
        self.assertEqual(1, ch3_cap_front.n_particles)
        self.assertEqual("_CH3", ch3_cap_front[0].name)
        self.assertEqual(1, len(ch3_cap_front.all_ports()))

    def test_dodecane(self):
        dodecane = AlkaneUA(12)
        self.assertEqual("_CH3", dodecane[0].name)
        self.assertEqual("_CH3", dodecane[-1].name)
        for i in range(1, 11):
            self.assertEqual("_CH2", dodecane[i].name)
        self.assertEqual(12, dodecane.n_particles)

    def test_n12_cap_front(self):
        n12_cap_front = AlkaneUA(n=12, cap_front=True, cap_end=False)
        self.assertEqual(12, n12_cap_front.n_particles)
        self.assertEqual("_CH3", n12_cap_front[0].name)
        self.assertEqual("_CH2", n12_cap_front[-1].name)
        self.assertEqual(1, len(n12_cap_front.all_ports()))
        self.assertIsInstance(n12_cap_front["down"], mb.Port)
        self.assertEqual(n12_cap_front[-1], n12_cap_front.all_ports()[0].anchor)

    def test_n12_cap_end(self):
        n12_cap_end = AlkaneUA(n=12, cap_front=False, cap_end=True)
        self.assertEqual(12, n12_cap_end.n_particles)
        self.assertEqual("_CH2", n12_cap_end[0].name)
        self.assertEqual("_CH3", n12_cap_end[-1].name)
        self.assertEqual(1, len(n12_cap_end.all_ports()))
        self.assertIsInstance(n12_cap_end["up"], mb.Port)
        self.assertEqual(n12_cap_end[0], n12_cap_end.all_ports()[0].anchor)

    def test_n12_no_caps(self):
        n12_no_cap = AlkaneUA(n=12, cap_front=False, cap_end=False)
        self.assertEqual(12, n12_no_cap.n_particles)
        self.assertEqual("_CH2", n12_no_cap[0].name)
        self.assertEqual("_CH2", n12_no_cap[-1].name)
        self.assertEqual(2, len(n12_no_cap.all_ports()))
        self.assertIsInstance(n12_no_cap["up"], mb.Port)
        self.assertIsInstance(n12_no_cap["down"], mb.Port)
        self.assertEqual(n12_no_cap[0], n12_no_cap.all_ports()[0].anchor)
        self.assertEqual(n12_no_cap[-1], n12_no_cap.all_ports()[1].anchor)


if __name__ == '__main__':
    unittest.main()
