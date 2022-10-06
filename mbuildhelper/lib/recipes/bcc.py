"""mBuild recipe for a BCC crystal."""
__all__ = ['BCC']

import mbuild as mb


class BCC(mb.Compound):
    """A BCC crystal of a single compound."""

    def __init__(self, lattice_spacing=None, compound_to_add=None, x=1, y=1, z=1, **kwargs):
        super(BCC, self).__init__(**kwargs)
        lattice_spacing = [lattice_spacing]*3
        lattice_points = {'A': [[0, 0, 0],
                                [0.5, 0.5, 0.5]]}
        angles = [90, 90, 90]

        fcc_lattice = mb.Lattice(lattice_spacing=lattice_spacing, angles=angles, lattice_points=lattice_points)
        self.add(fcc_lattice.populate(x=x, y=y, z=z, compound_dict={'A': compound_to_add}))


if __name__ == '__main__':
    m = BCC(lattice_spacing=0.285, compound_to_add=mb.Compound(name="Fe"), x=5, y=5, z=5)
    m.save("bcc_iron.mol2")
