"""A sp2 CH pseudoatom for united-atom force fields."""
__all__ = ['CHsp2UA']

import numpy as np

import mbuild as mb


# set the default bond length
DEFAULT_BOND_LENGTH = 0.133


class CHsp2UA(mb.Compound):
    """A united-atom CH pseudoatom."""

    def __init__(self, bond_length=DEFAULT_BOND_LENGTH, **kwargs):
        super(CHsp2UA, self).__init__(**kwargs)
        self.add(mb.Particle(name="_CH", pos=[0, 0, 0], mass=13.019, charge=0.0))

        self.add(mb.Port(anchor=self[0]), "up")
        self["up"].translate(np.array([0, bond_length/2, 0]))

        self.add(mb.Port(anchor=self[0]), "down")
        self["down"].translate(np.array([0, bond_length/2, 0]))
        self["down"].rotate(np.pi * 2 / 3, [0, 0, 1])


if __name__ == '__main__':
    m = CHsp2UA()
    m.save("ch_sp2_ua.mol2", overwrite=True)
