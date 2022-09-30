import mbuild as mb


class FCC(mb.Compound):

    def __init__(self, lattice_spacing=None, compound_to_add=None, x=1, y=1, z=1):
        super(FCC, self).__init__()
        lattice_points = {'A': [[0, 0, 0],
                                [0, 0.5, 0.5],
                                [0.5, 0.5, 0],
                                [0.5, 0, 0.5]]}
        angles = [90, 90, 90]

        fcc_lattice = mb.Lattice(lattice_spacing=lattice_spacing, angles=angles, lattice_points=lattice_points)
        self.add(fcc_lattice.populate(x=x, y=y, z=z, compound_dict={'A': compound_to_add}))
