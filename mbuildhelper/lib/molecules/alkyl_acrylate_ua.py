__all__ = ['AlkylAcrylateUA', 'ButylAcrylateUA',
           'LaurylAcrylateUA', 'DodecylAcrylateUA',
           'LaurylMethacrylateUA', 'DodecylMethacrylateUA']

import mbuild as mb

from mbuildhelper.lib.moieties import AcrylateUA, MethacrylateUA
import mbuildhelper.lib.recipes as recipes


class AlkylAcrylateUA(mb.Compound):

    def __init__(self, n=1, methyl=False, cap_front=True, cap_end=True, **kwargs):
        super(AlkylAcrylateUA, self).__init__(**kwargs)

        # add acrylate
        if methyl:
            self.add(MethacrylateUA(cap_front=cap_front, cap_end=cap_end), label="acrylate")
        else:
            self.add(AcrylateUA(cap_front=cap_front, cap_end=cap_end), label="acrylate")

        # add alkyl tail
        self.add(recipes.AlkaneUA(n=n, cap_front=False, cap_end=True), label="alkyl_tail")

        # create bond between alkyl tail and acrylate
        mb.force_overlap(
            move_this=self['alkyl_tail'],
            from_positions=self['alkyl_tail']['up'],
            to_positions=self['acrylate']['R']
        )

        # hoist ports to acrylate level
        if not cap_front:
            self.add(self['acrylate']['up'], "up", containment=False)
        if not cap_end:
            self.add(self['acrylate']['down'], "down", containment=False)


class ButylAcrylateUA(AlkylAcrylateUA):

    def __init__(self, cap_front=True, cap_end=True, **kwargs):
        super(ButylAcrylateUA, self).__init__(n=4, methyl=False, cap_front=cap_front, cap_end=cap_end, **kwargs)


class LaurylAcrylateUA(AlkylAcrylateUA):

    def __init__(self, cap_front=True, cap_end=True, **kwargs):
        super(LaurylAcrylateUA, self).__init__(n=12, methyl=False, cap_front=cap_front, cap_end=cap_end, **kwargs)


DodecylAcrylateUA = LaurylAcrylateUA


class LaurylMethacrylateUA(AlkylAcrylateUA):

    def __init__(self, cap_front=True, cap_end=True, **kwargs):
        super(LaurylMethacrylateUA, self).__init__(n=12, methyl=True, cap_front=cap_front, cap_end=cap_end, **kwargs)


DodecylMethacrylateUA = LaurylMethacrylateUA


if __name__ == '__main__':
    butyl_acrylate = AlkylAcrylateUA(n=4, cap_front=True, cap_end=True)
    butyl_acrylate.save("butyl_acrylate_ua.mol2", overwrite=True)
