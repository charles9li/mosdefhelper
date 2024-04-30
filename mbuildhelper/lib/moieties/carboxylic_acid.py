"""Carboxylic acid moiety."""
__all__ = ['CarboxylicAcid']

import mbuild as mb

from mbuildhelper.lib.moieties import Carbonyl, Hydroxyl


class CarboxylicAcid(mb.Compound):
    """A carboxylic acid group -C(=O)OH"""

    def __init__(self, **kwargs):
        super(CarboxylicAcid, self).__init__(**kwargs)

        # add hydroxyl and carbonyl groups
        self.add(Hydroxyl(), label="OH")
        self.add(Carbonyl(), label="C=O")

        # connect the groups
        mb.force_overlap(
            move_this=self['C=O'],
            from_positions=self['C=O']['R1'],
            to_positions=self['OH']['R']
        )

        # hoist port to carboxylic acid level
        self.add(self['C=O']['R2'], "R", containment=False)


if __name__ == '__main__':
    m = CarboxylicAcid()
    m.save('cooh.mol2', overwrite=True)
