"""A divalent, sp3 oxygen atom."""
__all__ = ['Osp3']

import numpy as np

import mbuild as mb


class Osp3(mb.Compound):
    """A divalent, sp3 oxygen.

    This atom can be used in both explicit-hydrogen and united-atom force fields.
    """

    def __init__(self, **kwargs):
        super(Osp3, self).__init__(**kwargs)
        self.add(mb.Particle(name="O", pos=[0, 0, 0], charge=0.0, element="O"))

        # angle to rotate ports
        theta = 0.5 * (180 - 115.0) * np.pi / 180

        # add ports
        self.add(mb.Port(anchor=self[0]), "up")
        self["up"].translate([0, -0.064, 0])
        self["up"].rotate(theta, around=[1, 0, 0])

        self.add(mb.Port(anchor=self[0]), "down")
        self["down"].translate([0, 0.064, 0])
        self["down"].rotate(np.pi, around=[0, 1, 0])
        self["down"].rotate(-theta, around=[1, 0, 0])


if __name__ == '__main__':
    m = Osp3()
    m.save("o_sp3.mol2")
