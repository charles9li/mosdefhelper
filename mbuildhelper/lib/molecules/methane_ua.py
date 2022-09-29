"""A methane pseudoatom for united-atom force fields."""
from mbuildhelper.lib.atoms import CH4UA


class MethaneUA(CH4UA):
    """A united-atom methane. Equivalent to the _CH4 pseudoatom."""
    pass


if __name__ == '__main__':
    m = MethaneUA()
    m.save("methane_ua.mol2", overwrite=True)
