"""mBuild recipe for a FCC crystal."""
__all__ = ['FCC']

import mbuild as mb


class FCC(mb.Compound):
    """A FCC crystal of a single compound."""

    def __init__(self, lattice_spacing=None, compound_to_add=None, x=1, y=1, z=1):
        super(FCC, self).__init__()
        lattice_spacing = [lattice_spacing]*3
        lattice_points = {'A': [[0, 0, 0],
                                [0, 0.5, 0.5],
                                [0.5, 0.5, 0],
                                [0.5, 0, 0.5]]}
        angles = [90, 90, 90]

        fcc_lattice = mb.Lattice(lattice_spacing=lattice_spacing, angles=angles, lattice_points=lattice_points)
        self.add(fcc_lattice.populate(x=x, y=y, z=z, compound_dict={'A': compound_to_add}))


if __name__ == '__main__':
    m = FCC(lattice_spacing=0.356, compound_to_add=mb.Compound(name="Fe"), x=5, y=5, z=5)
    m.save("fcc_iron.mol2")
