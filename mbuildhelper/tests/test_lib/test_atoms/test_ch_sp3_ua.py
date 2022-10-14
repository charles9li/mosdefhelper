import unittest

import mbuild as mb

from mbuildhelper.lib.atoms import CHsp3UA, CH3UA


class TestCHsp3UA(unittest.TestCase):

    def test_isobutane_ua(self):
        isobutane = mb.Compound()
        isobutane.add(CHsp3UA(), label="CH")
        for _ in range(3):
            isobutane.add(CH3UA(), label="CH3[$]")
        for i, port in enumerate(['up', 'down', 'left']):
            mb.force_overlap(
                move_this=isobutane[f'CH3[{i}]'],
                from_positions=isobutane[f'CH3[{i}]']['up'],
                to_positions=isobutane['CH'][port]
            )
        isobutane.save('isobutane_ua.mol2', overwrite=True)


if __name__ == '__main__':
    unittest.main()
