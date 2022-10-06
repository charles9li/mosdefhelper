"""A CH4 (methane) pseudoatom for united-atom force fields."""
__all__ = ['CH4UA']

import mbuild as mb


class CH4UA(mb.Compound):
    """A united-atom methane."""

    def __init__(self, **kwargs):
        super(CH4UA, self).__init__(**kwargs)
        self.add(mb.Particle(name='_CH4', pos=[0, 0, 0], mass=16.043))


if __name__ == '__main__':
    m = CH4UA()
    m.save("ch4_ua.mol2", overwrite=True)
