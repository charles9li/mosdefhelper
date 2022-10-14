import unittest

import mbuild as mb

from mbuildhelper.lib.atoms import CH3UA
from mbuildhelper.lib.moieties import Carbonyl


class TestCarbonyl(unittest.TestCase):

    def test_acetone_ua(self):
        acetone = mb.Compound()
        acetone.add(Carbonyl(), label="C=O")
        acetone.add(CH3UA(), label="CH3[$]")
        acetone.add(CH3UA(), label="CH3[$]")
        print(acetone.all_ports())
        mb.force_overlap(
            move_this=acetone['CH3[0]'],
            from_positions=acetone['CH3[0]']['up'],
            to_positions=acetone['C=O']['R1']
        )
        mb.force_overlap(
            move_this=acetone['CH3[1]'],
            from_positions=acetone['CH3[1]']['up'],
            to_positions=acetone['C=O']['R2']
        )
        acetone.save("acetone_ua.mol2", overwrite=True)


if __name__ == '__main__':
    unittest.main()
