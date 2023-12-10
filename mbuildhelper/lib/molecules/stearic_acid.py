"""A stearic acid molecule for united-atom force fields."""
__all__ = ['StearicAcid']

import mbuild as mb
from mbuild.lib.recipes import Alkane

from mbuildhelper.lib.moieties import CarboxylicAcid


class StearicAcid(mb.Compound):
    """A united-atom stearic acid molecule."""

    def __init__(self, **kwargs):
        super(StearicAcid, self).__init__(**kwargs)

        self.add(mb.load("OC(=O)CCCCCCCCCCCCCCCCC", smiles=True))

        # # add carboxylic acid group
        # self.add(CarboxylicAcid(), label="COOH")
        #
        # # add C17 tail
        # self.add(Alkane(n=17, cap_front=False, cap_end=True), label="C17")
        #
        # # create bond between C17 tail and COOH
        # mb.force_overlap(
        #     move_this=self['C17'],
        #     from_positions=self['C17']['up'],
        #     to_positions=self['COOH']['R']
        # )


if __name__ == '__main__':
    stearic_acid = StearicAcid()
    stearic_acid.save('stearic_acid.mol2', overwrite=True)
