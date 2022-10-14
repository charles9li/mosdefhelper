"""Ester moiety."""
__all__ = ['Ester']

import mbuild as mb

from mbuildhelper.lib.atoms import Osp3
from mbuildhelper.lib.moieties import Carbonyl


class Ester(mb.Compound):
    """A ester group -C(=O)O-."""

    def __init__(self, **kwargs):
        super(Ester, self).__init__(**kwargs)

        self.add(Carbonyl(), label="C=O")   # carbonyl
        self.add(Osp3(), label="O[$]")      # ether oxygen
        mb.force_overlap(
            move_this=self['O[0]'],
            from_positions=self['O[0]']['up'],
            to_positions=self['C=O']['R2']
        )

        # hoist ports to ester level
        self.add(self['C=O']['R1'], "R1", containment=False)
        self.add(self['O[0]']['down'], "R2", containment=False)


if __name__ == '__main__':
    m = Ester()
    m.save('ester.mol2', overwrite=True)
