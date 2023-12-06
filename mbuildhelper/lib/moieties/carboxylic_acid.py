"""Carboxylic acid moiety."""
__all__ = ['CarboxylicAcid']

import mbuild as mb

from mbuildhelper.lib.atoms import Osp3
from mbuildhelper.lib.moieties import Carbonyl


class _H(mb.Compound):
    """Acidic hydrogen in carboxylic acid group."""

    def __init__(self, **kwargs):
        super(_H, self).__init__(**kwargs)
        self.add(mb.Particle(name="H", pos=[0, 0, 0], charge=0.0, element="H"))

        # add port
        self.add(mb.Port(anchor=self[0]), label="up")
        self["up"].translate([0, 0.033, 0])


class CarboxylicAcid(mb.Compound):
    """A carboxylic acid group -C(=O)OH"""

    def __init__(self, **kwargs):
        super(CarboxylicAcid, self).__init__(**kwargs)

        # add hydroxyl group
        self.add(Osp3(), label="O[$]")
        self.add(_H(), label="H[$]")
        mb.force_overlap(
            move_this=self['O[0]'],
            from_positions=self['O[0]']['up'],
            to_positions=self['H[0]']['up']
        )

        # add carbonyl group
        self.add(Carbonyl(), label="C=O")
        mb.force_overlap(
            move_this=self['C=O'],
            from_positions=self['C=O']['R1'],
            to_positions=self['O[0]']['down']
        )

        # hoist port to carboxylic acid level
        self.add(self['C=O']['R2'], "R", containment=False)


if __name__ == '__main__':
    m = CarboxylicAcid()
    m.save('cooh.mol2', overwrite=True)
