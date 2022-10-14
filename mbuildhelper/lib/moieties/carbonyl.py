"""Carbonyl moiety."""
__all__ = ['Carbonyl']

import mbuild as mb

from mbuildhelper.lib.atoms import Csp2, Osp2


class Carbonyl(mb.Compound):
    """A carbonyl group RR'C(=O)"""

    def __init__(self, **kwargs):
        super(Carbonyl, self).__init__(**kwargs)

        # add carbon and oxygen
        self.add(Csp2(), label="C[$]")
        self.add(Osp2(), label="O[$]")

        # create double bond
        mb.force_overlap(
            move_this=self['O[0]'],
            from_positions=self['O[0]']['up'],
            to_positions=self['C[0]']['up']
        )

        # hoist ports up to carbonyl level
        self.add(self['C[0]']['down'], "R1", containment=False)
        self.add(self['C[0]']['left'], "R2", containment=False)


if __name__ == '__main__':
    m = Carbonyl()
    m.save('carbonyl.mol2', overwrite=True)
