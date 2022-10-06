"""A methane pseudoatom for united-atom force fields."""
__all__ = ['MethaneUA']

from mbuildhelper.lib.atoms import CH4UA


class MethaneUA(CH4UA):
    """A united-atom methane. Equivalent to the _CH4 pseudoatom."""

    def __init__(self, **kwargs):
        super(MethaneUA, self).__init__(**kwargs)


if __name__ == '__main__':
    m = MethaneUA()
    m.save("methane_ua.mol2", overwrite=True)
