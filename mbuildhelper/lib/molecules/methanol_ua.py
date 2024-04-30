"""A methanol molecule for united-atom force fields."""
__all__ = ['MethanolUA']

import mbuild as mb

from mbuildhelper.lib.atoms import CH3UA
from mbuildhelper.lib.moieties import Hydroxyl


class MethanolUA(mb.Compound):
    """A united-atom methanol molecule."""

    def __init__(self, **kwargs):
        super(MethanolUA, self).__init__(**kwargs)

        # add hydroxyl group and CH3 pseudoatom
        self.add(Hydroxyl(), label="OH")
        self.add(CH3UA(), label="CH3")

        # create bond between hydroxyl and OH
        mb.force_overlap(
            move_this=self['OH'],
            from_positions=self['OH']['R'],
            to_positions=self['CH3']['up']
        )


if __name__ == '__main__':
    methanol = MethanolUA()
    methanol.save('methanol_ua.mol2', overwrite=True)
