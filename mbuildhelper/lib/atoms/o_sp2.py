"""A sp2 oxygen atom."""
__all__ = ['Osp2']

import numpy as np

import mbuild as mb


class Osp2(mb.Compound):
    """A sp2 oxygen.

    This atom can be used in both explicit-hydrogen and united-atom force fields.
    """

    def __init__(self, **kwargs):
        super(Osp2, self).__init__(**kwargs)
        self.add(mb.Particle(name="O", pos=[0, 0, 0], charge=0.0, element="O"))

        # add ports
        self.add(mb.Port(anchor=self[0]), "up")
        self["up"].translate([0, 0.1200/2, 0])


if __name__ == '__main__':
    m = Osp2()
    m.save("o_sp2.mol2")
