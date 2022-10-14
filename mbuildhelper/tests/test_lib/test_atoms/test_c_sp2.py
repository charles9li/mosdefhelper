import unittest

import mbuild as mb

from mbuildhelper.lib.atoms import Csp2, CH3UA


class TestCsp2(unittest.TestCase):

    def test_angles(self):
        c = mb.Compound()
        c.add(Csp2())
        for port in ["up", "down", "left"]:
            ch3 = CH3UA()
            c.add(ch3)
            mb.force_overlap(
                move_this=ch3,
                from_positions=ch3["up"],
                to_positions=c["Csp2[0]"][port]
            )
        c.save("test_angles.mol2", overwrite=True)
