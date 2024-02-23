import os
from pkg_resources import resource_filename

import numpy as np

import mbuild as mb
from mbuild.lattice import load_cif
from mbuild.lib.recipes import TiledCompound


class Hematite0001SurfaceTiled(mb.Compound):

    def __init__(self, n_tiles=None, lengths=None, add_linear_bonds=False, **kwargs):
        super(Hematite0001SurfaceTiled, self).__init__(**kwargs)

        # load in data from cif file
        cif_filepath = os.path.join(resource_filename("mbuildhelper", "lib"), "surfaces", "fe2o3.cif")
        fe2o3_lattice = load_cif(cif_filepath)

        # create atoms
        fe = mb.Compound(name='Fe', element='Fe')
        o = mb.Compound(name='O', element='O')

        # populate lattice and create rectilinear unit cell
        compound_dict = {'O': o, 'Fe': fe}
        fe2o3_unit_cell = fe2o3_lattice.populate(compound_dict=compound_dict, x=1, y=2, z=1)

        # modify unit cell such that it is rectangular in xy plane
        new_box = mb.Box(lengths=np.diagonal(fe2o3_unit_cell.box.vectors))
        fe2o3_unit_cell.box = new_box

        # shift atom positions into new box
        fe2o3_unit_cell.xyz -= new_box.lengths * np.round_(fe2o3_unit_cell.xyz / new_box.lengths - 0.5)
        fe2o3_unit_cell.periodicity = (True, True, True)

        # check that either n_tiles or lengths are specified
        if n_tiles is None and lengths is None:
            raise ValueError("n_tiles or lengths must be specified")
        if n_tiles is not None and lengths is not None:
            raise ValueError("only one of n_tiles or lengths can be specified")

        # if lengths is specified, determine the number of tiles using the unit cell dimensions
        if lengths is not None:
            _lattice_spacing = np.array(fe2o3_unit_cell.box.lengths) / 10.
            n_tiles = [int(_l / _lattice_spacing[i]) + 1 for i, _l in enumerate(lengths)]

        # tile rectilinear unit cell
        fe2o3_tiled = TiledCompound(fe2o3_unit_cell, n_tiles=n_tiles)

        # add bonds (needed for foyer parameterization)
        prev_particle = None
        if add_linear_bonds:
            for particle in fe2o3_tiled.particles():
                if prev_particle is not None:
                    fe2o3_tiled.add_bond((prev_particle, particle))
                prev_particle = particle
        else:
            # get particles
            particles = list(fe2o3_tiled.particles())

            # compute square pair distances
            xyz = fe2o3_tiled.xyz
            rij = xyz - xyz[:, None, :]
            lengths = np.array(new_box.lengths)
            rij = rij - np.array([1, 1, 0]) * lengths * np.round(rij / lengths)
            rij_sqd = np.sum(rij**2, axis=2)

            # determine which pairs are within 3 angstrom cutoff
            rij_sqd[np.tril_indices_from(rij_sqd)] = np.inf # fill lower triangle of square dist matrix with inf to prevent double counting pairs
            pairs = np.transpose(np.where(rij_sqd < .3**2))

            # add bonds for all pairs within cutoff
            for p in pairs:
                fe2o3_tiled.add_bond((particles[p[0]], particles[p[1]]))

        # add lattice
        self.add(fe2o3_tiled)


if __name__ == '__main__':
    hematite = Hematite0001SurfaceTiled(lengths=[3, 3, 2], add_linear_bonds=True)
    hematite.save("hematite_0001_tiled.pdb", overwrite=True)
