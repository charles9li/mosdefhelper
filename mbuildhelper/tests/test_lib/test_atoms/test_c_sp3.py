import unittest

import mbuild as mb

from mbuildhelper.lib.atoms import Csp3, CH3UA


class TestCsp3UA(unittest.TestCase):

    def test_neopentaen_ua(self):
        neopentane = mb.Compound()
        neopentane.add(Csp3(), label="C")
        for _ in range(4):
            neopentane.add(CH3UA(), label="CH3[$]")
        for i, port in enumerate(['up', 'down', 'left', 'right']):
            mb.force_overlap(
                move_this=neopentane[f'CH3[{i}]'],
                from_positions=neopentane[f'CH3[{i}]']['up'],
                to_positions=neopentane['C'][port]
            )
        neopentane.save('neopentane_ua.mol2', overwrite=True)


if __name__ == '__main__':
    unittest.main()
