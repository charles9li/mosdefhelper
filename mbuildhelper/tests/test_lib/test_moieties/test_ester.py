import unittest

import mbuild as mb

from mbuildhelper.lib.atoms import CH3UA
from mbuildhelper.lib.moieties import Ester


class TestEster(unittest.TestCase):

    def test_methyl_acetate_ua(self):
        methyl_acetate = mb.Compound()
        methyl_acetate.add(Ester(), label="OC=O")
        methyl_acetate.add(CH3UA(), label="CH3[$]")
        methyl_acetate.add(CH3UA(), label="CH3[$]")
        mb.force_overlap(
            move_this=methyl_acetate['CH3[0]'],
            from_positions=methyl_acetate['CH3[0]']['up'],
            to_positions=methyl_acetate['OC=O']['R1']
        )
        mb.force_overlap(
            move_this=methyl_acetate['CH3[1]'],
            from_positions=methyl_acetate['CH3[1]']['up'],
            to_positions=methyl_acetate['OC=O']['R2']
        )
        methyl_acetate.save("methyl_acetate_ua.mol2", overwrite=True)


if __name__ == '__main__':
    unittest.main()
