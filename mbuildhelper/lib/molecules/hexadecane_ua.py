"""A dodecane molecule for united-atom force fields."""
__all__ = ['HexadecaneUA']

import mbuild as mb

import mbuildhelper.lib.recipes as mbh_lib_recipes


class HexadecaneUA(mb.Compound):
    """A dodecane molecule for united-atom force fields.

    Uses the recipe for united-atom alkanes.
    """

    def __init__(self, **kwargs):
        super(HexadecaneUA, self).__init__(**kwargs)
        self.add(mbh_lib_recipes.AlkaneUA(n=16, cap_front=True, cap_end=True, **kwargs))


if __name__ == '__main__':
    m = HexadecaneUA()
    m.save("hexadecane_ua.mol2", overwrite=True)
