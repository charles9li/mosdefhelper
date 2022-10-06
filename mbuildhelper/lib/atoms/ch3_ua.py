"""A CH3 pseudoatom for united-atom force fields."""
__all__ = ['CH3UA']

import mbuild as mb


class CH3UA(mb.Compound):
    """A united-atom methyl group."""

    def __init__(self, **kwargs):
        super(CH3UA, self).__init__(**kwargs)
        self.add(mb.Particle(name='_CH3', pos=[0, 0, 0], mass=15.035))
        self.add(mb.Port(anchor=self[0]), "up")
        self["up"].translate([0, -0.154/2, 0])


if __name__ == '__main__':
    m = CH3UA()
    m.save("ch3_ua.mol2", overwrite=True)
