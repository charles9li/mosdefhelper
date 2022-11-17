"""A CH2 sp2 pseudoatom for united-atom force fields."""
__all__ = ['CH2sp2UA']

import mbuild as mb


class CH2sp2UA(mb.Compound):
    """A CH2 sp2 pseudoatom for united-atom force fields."""

    def __init__(self, **kwargs):
        super(CH2sp2UA, self).__init__(**kwargs)
        self.add(mb.Particle(name='_CH2', pos=[0, 0, 0], mass=14.027, charge=0.0))
        self.add(mb.Port(anchor=self[0]), "up")
        self["up"].translate([0, -0.133/2, 0])


if __name__ == '__main__':
    m = CH2sp2UA()
    m.save("ch2_sp2_ua.mol2", overwrite=True)
