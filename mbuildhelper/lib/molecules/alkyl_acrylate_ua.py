__all__ = ['AlkylAcrylate', 'ButylAcrylate',
           'LaurylAcrylate', 'DodecylAcrylate',
           'LaurylMethacrylate', 'DodecylMethacrylate']

import mbuild as mb

from mbuildhelper.lib.moieties import AcrylateUA, MethacrylateUA
from mbuildhelper.lib.recipes import AlkaneUA


class AlkylAcrylate(mb.Compound):

    def __init__(self, n=1, methyl=False, cap_front=True, cap_end=True, **kwargs):
        super(AlkylAcrylate, self).__init__(**kwargs)

        # add acrylate
        if methyl:
            self.add(MethacrylateUA(cap_front=cap_front, cap_end=cap_end), label="acrylate")
        else:
            self.add(AcrylateUA(cap_front=cap_front, cap_end=cap_end), label="acrylate")

        # add alkyl tail
        self.add(AlkaneUA(n=n, cap_front=False, cap_end=True), label="alkyl_tail")

        # create bond between alkyl tail and acrylate
        mb.force_overlap(
            move_this=self['alkyl_tail'],
            from_positions=self['alkyl_tail']['up'],
            to_positions=self['acrylate']['R']
        )


class ButylAcrylate(AlkylAcrylate):

    def __init__(self, cap_front=True, cap_end=True, **kwargs):
        super(ButylAcrylate, self).__init__(n=4, methyl=False, cap_front=cap_front, cap_end=cap_end, **kwargs)


class LaurylAcrylate(AlkylAcrylate):

    def __init__(self, cap_front=True, cap_end=True, **kwargs):
        super(LaurylAcrylate, self).__init__(n=12, methyl=False, cap_front=cap_front, cap_end=cap_end, **kwargs)


DodecylAcrylate = LaurylAcrylate


class LaurylMethacrylate(AlkylAcrylate):

    def __init__(self, cap_front=True, cap_end=True, **kwargs):
        super(LaurylMethacrylate, self).__init__(n=12, methyl=True, cap_front=cap_front, cap_end=cap_end, **kwargs)


DodecylMethacrylate = LaurylMethacrylate


if __name__ == '__main__':
    butyl_acrylate = AlkylAcrylate(n=4, cap_front=True, cap_end=True)
    butyl_acrylate.save("butyl_acrylate_ua.mol2", overwrite=True)
