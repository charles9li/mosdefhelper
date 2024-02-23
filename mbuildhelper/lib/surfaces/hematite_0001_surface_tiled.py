"""Hematite (0001) surface."""
__all__ = ['Hematite0001SurfaceTiled']

import numpy as np

import mbuild as mb

from mbuildhelper.lib.recipes import SimpleTiledCompound
from mbuildhelper.lib.surfaces import Hematite0001RectilinearUnitCell


class Hematite0001SurfaceTiled(mb.Compound):
    """Hematite (0001) surface.

    The surface is created using a custom SimpleTiledCompound recipe instead of
    the populate method of the Lattice class.
    """

    def __init__(self, n_tiles=None, lengths=None, add_bonds=False, **kwargs):
        super(Hematite0001SurfaceTiled, self).__init__(**kwargs)

        # check that either n_tiles or lengths are specified
        if n_tiles is None and lengths is None:
            raise ValueError("n_tiles or lengths must be specified")
        if n_tiles is not None and lengths is not None:
            raise ValueError("only one of n_tiles or lengths can be specified")

        # load in unit cell
        fe2o3_unit_cell = Hematite0001RectilinearUnitCell(add_bonds=add_bonds)
        fe2o3_unit_cell.periodicity = (True, True, True)

        # if lengths is specified, determine the number of tiles using the unit cell dimensions
        if lengths is not None:
            _lattice_spacing = np.array(fe2o3_unit_cell.box.lengths)
            n_tiles = [int(_l / _lattice_spacing[i]) + 1 for i, _l in enumerate(lengths)]

        # tile rectilinear unit cell
        fe2o3_tiled = SimpleTiledCompound(fe2o3_unit_cell, n_tiles=n_tiles)

        # add lattice
        self.add(fe2o3_tiled)


if __name__ == '__main__':
    hematite = Hematite0001SurfaceTiled(lengths=[2, 2, 1], add_bonds=True)
    hematite.save("hematite_0001_tiled.pdb", overwrite=True)
