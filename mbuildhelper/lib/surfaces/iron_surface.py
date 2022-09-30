import mbuild as mb

from mbuildhelper.lib.recipes import BCC, FCC


_CRYSTAL_STRUCTURE_TO_LATTICE_DICT = {'bcc': BCC,
                                      'fcc': FCC}


def _compute_fe_lattice_constant(crystal_structure, temperature):
    if crystal_structure == 'bcc':
        return 4.393e-4 * temperature + 0.28473
    elif crystal_structure == 'fcc':
        return 5.84340565e-06 * temperature + 3.55526286e-01
    else:
        raise ValueError("invalid crystal structure")


class IronSurface(mb.Compound):
    """Creates an iron surface.

    Parameters
    ----------
    crystal_structure : str
        BCC or FCC
    temperature : float
        Temperature of the crystal. Used for calculating lattice constants.
    n_tiles : array-like, shape=(3,), dtype=int
        Number of times to replicate title in the x, y, and z-directions.
    lengths : array-like, shape=(3,), dtype=float
        Dimensions of slab in nanometers in each direction.
    """

    def __init__(self, crystal_structure=None, temperature=None, n_tiles=None, lengths=None):
        super(IronSurface, self).__init__()

        # check that crystal structure is valid
        crystal_structure = crystal_structure.lower()
        if crystal_structure not in ['bcc', 'fcc']:
            raise ValueError("invalid crystal structure")

        # get lattice constant
        lattice_spacing = _compute_fe_lattice_constant(crystal_structure, temperature)

        # check that either n_tiles or lengths are specified
        if n_tiles is None and lengths is None:
            raise ValueError("n_tiles or lengths must be specified")
        if n_tiles is not None and lengths is not None:
            raise ValueError("only one of n_tiles or lengths can be specified")
        if lengths is not None:
            n_tiles = [int(_l / lattice_spacing) + 1 for _l in lengths]

        # create lattice
        lattice = _CRYSTAL_STRUCTURE_TO_LATTICE_DICT[crystal_structure](lattice_spacing=lattice_spacing,
                                                                        compound_to_add=mb.Compound(name="Fe"),
                                                                        x=n_tiles[0], y=n_tiles[1], z=n_tiles[2])
        self.add(lattice)


if __name__ == '__main__':
    surface = IronSurface(crystal_structure='fcc', temperature=313.15, lengths=[5, 5, 2])
    print(surface)
