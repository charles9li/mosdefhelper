"""Hydroxyl moiety."""
__all__ = ['Hydroxyl']

import mbuild as mb

from mbuildhelper.lib.atoms import Osp3


class _H(mb.Compound):
    """Acidic hydrogen in hydroxyl group."""

    def __init__(self, **kwargs):
        super(_H, self).__init__(**kwargs)
        self.add(mb.Particle(name="H", pos=[0, 0, 0], charge=0.0, element="H"))

        # add port
        self.add(mb.Port(anchor=self[0]), label="up")
        self["up"].translate([0, 0.033, 0])


class Hydroxyl(mb.Compound):
    """A hydroxyl group OH"""

    def __init__(self, **kwargs):
        super(Hydroxyl, self).__init__(**kwargs)

        # add hydroxyl group
        self.add(Osp3(), label="O[$]")
        self.add(_H(), label="H[$]")
        mb.force_overlap(
            move_this=self['O[0]'],
            from_positions=self['O[0]']['up'],
            to_positions=self['H[0]']['up']
        )

        # hoist port to hydroxyl level
        self.add(self['O[0]']['down'], "R", containment=False)


if __name__ == '__main__':
    m = Hydroxyl()
    m.save('oh.mol2', overwrite=True)
