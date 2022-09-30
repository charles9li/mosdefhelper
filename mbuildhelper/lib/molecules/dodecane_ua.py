"""A dodecane molecule for united-atom force fields."""
from mbuildhelper.lib.recipes import AlkaneUA


class DodecaneUA(AlkaneUA):
    """A dodecane molecule for united-atom force fields.

    Uses the recipe for united-atom alkanes.
    """

    def __init__(self):
        super(DodecaneUA, self).__init__(n=12, cap_front=True, cap_end=True)


if __name__ == '__main__':
    m = DodecaneUA()
    m.save("dodecane_ua.mol2", overwrite=True)
