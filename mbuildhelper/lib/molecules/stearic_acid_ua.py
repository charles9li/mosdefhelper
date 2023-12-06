"""A stearic acid molecule for united-atom force fields."""
__all__ = ['StearicAcidUA']

import mbuild as mb

from mbuildhelper.lib.moieties import CarboxylicAcid
from mbuildhelper.lib.recipes import AlkaneUA


class StearicAcidUA(mb.Compound):
    """A united-atom stearic acid molecule."""

    def __init__(self, **kwargs):
        super(StearicAcidUA, self).__init__(**kwargs)

        # add carboxylic acid group
        self.add(CarboxylicAcid(), label="COOH")

        # add C17 tail
        self.add(AlkaneUA(n=17, cap_front=False, cap_end=True), label="C17")

        # create bond between C17 tail and COOH
        mb.force_overlap(
            move_this=self['C17'],
            from_positions=self['C17']['up'],
            to_positions=self['COOH']['R']
        )


if __name__ == '__main__':
    stearic_acid = StearicAcidUA()
    stearic_acid.save('stearic_acid_ua.mol2', overwrite=True)
