"""A dodecane molecule for united-atom force fields."""
__all__ = ['DodecaneUA']

import mbuild as mb

import mbuildhelper.lib.recipes as mbh_lib_recipes


class DodecaneUA(mb.Compound):
    """A dodecane molecule for united-atom force fields.

    Uses the recipe for united-atom alkanes.
    """

    def __init__(self, **kwargs):
        super(DodecaneUA, self).__init__(**kwargs)
        self.add(mbh_lib_recipes.AlkaneUA(n=12, cap_front=True, cap_end=True, **kwargs))


if __name__ == '__main__':
    m = DodecaneUA()
    m.save("dodecane_ua.mol2", overwrite=True)
