"""United-atom phenyl moiety."""
__all__ = ['PhenylUA6Site', 'PhenylUA9Site']

import numpy as np

import mbuild as mb


SP2_BOND_LENGTH = 0.14
SP3_BOND_LENGTH = 0.154


class PhenylUA6Site(mb.Compound):
    """A 6-atom version of a phenyl group"""

    def __init__(self, cap=True, **kwargs):
        super(PhenylUA6Site, self).__init__(**kwargs)

        # add first carbon
        pos = [0, SP2_BOND_LENGTH, 0]
        if cap:
            self.add(mb.Particle(name="_CH", pos=pos, charge=0.0, mass=13.019))
        else:
            self.add(mb.Particle(name="_C", pos=pos, charge=0.0, element="C"))
            # add port
            self.add(mb.Port(anchor=self[0]), "R")
            self["R"].translate(np.array([0, SP3_BOND_LENGTH/2, 0]))

        # add other carbons in the ring
        for i in range(1, 6):
            self.add(mb.Particle(name="_CH", pos=pos, charge=0.0, mass=13.019))
            self[i].rotate(np.pi * i / 3, [0, 0, 1])

        # add bonds between the carbons in the ring
        for i in range(5):
            self.add_bond([self[i], self[i+1]])
        self.add_bond([self[5], self[0]])


# distance between center and pi quadrupole sites
MCENTER_MPI_BOND_LENGTH = 0.0785


class PhenylUA9Site(PhenylUA6Site):
    """A 9-atom version of a phenyl group"""

    def __init__(self, cap=True, **kwargs):
        # call __init__ of 6-site UA phenyl
        super(PhenylUA9Site, self).__init__(cap=cap, **kwargs)

        # add virtual sites
        m_center = mb.Particle(
            name="_MC",
            pos=[0, 0, 0],
            mass=0.0,
            charge=0.242,
            element=None
        )
        m_pi_top = mb.Particle(
            name="_MP",
            pos=[0, 0, MCENTER_MPI_BOND_LENGTH],
            mass=0.0,
            charge=-0.121,
            element=None
        )
        m_pi_bottom = mb.Particle(
            name="_MP",
            pos=[0, 0, -MCENTER_MPI_BOND_LENGTH],
            mass=0.0,
            charge=-0.121,
            element=None
        )
        self.add([m_center, m_pi_top, m_pi_bottom])


if __name__ == '__main__':
    phenyl_ua_6site = PhenylUA6Site()
    phenyl_ua_6site.save("phenyl_ua_6site.mol2", overwrite=True)
    phenyl_ua_9site = PhenylUA9Site()
    phenyl_ua_9site.save("phenyl_ua_9site.mol2", overwrite=True)
