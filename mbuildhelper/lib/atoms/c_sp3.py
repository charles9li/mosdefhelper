"""A sp3 C pseudoatom for united-atom force fields."""
__all__ = ['Csp3']

import numpy as np

import mbuild as mb


class Csp3(mb.Compound):
    """A united-atom methylene linker."""

    def __init__(self, **kwargs):
        super(Csp3, self).__init__(**kwargs)
        self.add(mb.Particle(name='C', pos=[0, 0, 0], element='C'))

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

        self.add(mb.Port(anchor=self[0], orientation=[-1, 0, 0], separation=0.154/2), "left")
        self["left"].rotate(theta, around=[0, 1, 0])

        self.add(mb.Port(anchor=self[0], orientation=[1, 0, 0], separation=0.154/2), "right")
        self["right"].rotate(-theta, around=[0, 1, 0])


if __name__ == '__main__':
    m = Csp3()
    m.save("c_sp3_ua.mol2")
