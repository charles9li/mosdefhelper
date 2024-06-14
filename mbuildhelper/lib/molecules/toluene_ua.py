"""A toluene molecule for united-atom force fields."""
__all__ = ['TolueneUA7Site', 'TolueneUA10Site']

import mbuild as mb

from mbuildhelper.lib.atoms import CH3UA
from mbuildhelper.lib.moieties import PhenylUA6Site, PhenylUA9Site


class TolueneUA7Site(mb.Compound):
    """A 6-site united-atom toluene molecule."""

    def __init__(self, **kwargs):
        super(TolueneUA7Site, self).__init__(**kwargs)

        # add phenyl moiety
        self.add(PhenylUA6Site(cap=False), label="phenyl")

        # add methyl group
        self.add(CH3UA(), label="methyl")

        # create bond between phenyl and methyl groups
        mb.force_overlap(
            move_this=self["methyl"],
            from_positions=self["methyl"]["up"],
            to_positions=self["phenyl"]["R"]
        )


class TolueneUA10Site(mb.Compound):
    """A 9-site united-atom toluene molecule."""

    def __init__(self, **kwargs):
        super(TolueneUA10Site, self).__init__(**kwargs)

        # add phenyl moiety
        self.add(PhenylUA9Site(cap=False), label="phenyl")

        # add methyl group
        self.add(CH3UA(), label="methyl")

        # create bond between phenyl and methyl groups
        mb.force_overlap(
            move_this=self["methyl"],
            from_positions=self["methyl"]["up"],
            to_positions=self["phenyl"]["R"]
        )


if __name__ == '__main__':
    toluene_ua_7site = TolueneUA10Site()
    toluene_ua_7site.save("toluene_ua_7site.mol2", overwrite=True)
    toluene_ua_10site = TolueneUA10Site()
    toluene_ua_10site.save("toluene_ua_10site.mol2", overwrite=True)
