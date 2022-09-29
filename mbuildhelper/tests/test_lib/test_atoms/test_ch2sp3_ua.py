import unittest

from mbuildhelper.lib.atoms import CH2sp3UA


class TestCH2sp3UA(unittest.TestCase):

    def test_atom(self):
        m = CH2sp3UA()
        self.assertEqual(2, len(m.all_ports()))
        self.assertEqual(2, len(m.available_ports()))
