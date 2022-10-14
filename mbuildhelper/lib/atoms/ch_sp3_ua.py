"""A sp3 CH pseudoatom for united-atom force fields."""
__all__ = ['CHsp3UA']

import numpy as np

import mbuild as mb


class CHsp3UA(mb.Compound):
    """A united-atom methylene linker."""

    def __init__(self, **kwargs):
        super(CHsp3UA, self).__init__(**kwargs)
        self.add(mb.Particle(name='_CH', pos=[0, 0, 0], mass=13.019))

        # angle to rotate ports
        theta = 0.5 * (180 - 109.5) * np.pi / 180

        # add ports
        self.add(mb.Port(anchor=self[0]), "up")
        self['up'].translate([0, -0.154/2, 0])
        self["up"].rotate(theta, around=[1, 0, 0])

        self.add(mb.Port(anchor=self[0]), "down")
        self['down'].translate([0, 0.154/2, 0])
        self["down"].rotate(np.pi, around=[0, 1, 0])
        self["down"].rotate(-theta, around=[1, 0, 0])

        self.add(mb.Port(anchor=self[0], orientation=[-1, 0, 0], separation=.154/2), "left")
        self["left"].rotate(theta, around=[0, 1, 0])


if __name__ == '__main__':
    m = CHsp3UA()
    m.save("ch_sp3_ua.mol2")
