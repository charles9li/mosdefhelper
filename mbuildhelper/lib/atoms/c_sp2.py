"""A sp2 C atom."""
__all__ = ['Csp2']

import numpy as np

import mbuild as mb


# set the default bond length
DEFAULT_BOND_LENGTH = 0.133


class Csp2(mb.Compound):
    """A sp2 C atom.
    
    This atom can be used in both explicit-hydrogen and united-atom force fields.
    """
    
    def __init__(self, bond_length=DEFAULT_BOND_LENGTH, **kwargs):
        super(Csp2, self).__init__(**kwargs)
        self.add(mb.Particle(name="C", pos=[0, 0, 0], charge=0.0, element="C"))

        self.add(mb.Port(anchor=self[0]), "up")
        self["up"].translate(np.array([0, bond_length/2, 0]))

        self.add(mb.Port(anchor=self[0]), "down")
        self["down"].translate(np.array([0, bond_length/2, 0]))
        self["down"].rotate(np.pi * 2 / 3, [0, 0, 1])

        self.add(mb.Port(anchor=self[0]), "left")
        self["left"].translate(np.array([0, bond_length/2, 0]))
        self["left"].rotate(-np.pi * 2 / 3, [0, 0, 1])


if __name__ == '__main__':
    m = Csp2()
    m.save("c_sp2.mol2", overwrite=True)
