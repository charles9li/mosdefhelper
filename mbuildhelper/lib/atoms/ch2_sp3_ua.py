"""A sp3 CH2 pseudoatom for united-atom force fields."""
__all__ = ['CH2sp3UA']

import numpy as np

import mbuild as mb


class CH2sp3UA(mb.Compound):
    """A united-atom methylene linker."""

    def __init__(self, **kwargs):
        super(CH2sp3UA, self).__init__(**kwargs)
        self.add(mb.Particle(name='_CH2', pos=[0, 0, 0], mass=14.027, charge=0.0))

        # angle to rotate ports
        theta = 0.5 * (180 - 109.5) * np.pi / 180

        # add ports
        self.add(mb.Port(anchor=self[0]), "up")
        self["up"].translate([0, -0.154/2, 0])
        self["up"].rotate(theta, around=[1, 0, 0])

        self.add(mb.Port(anchor=self[0]), "down")
        self["down"].translate([0, 0.154/2, 0])
        self["down"].rotate(np.pi, around=[0, 1, 0])
        self["down"].rotate(-theta, around=[1, 0, 0])


if __name__ == '__main__':
    m = CH2sp3UA()
    m.save("ch2_sp3_ua.mol2", overwrite=True)
