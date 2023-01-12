import os
from pkg_resources import resource_filename

import numpy as np

import mbuild as mb
from mbuild.lattice import load_cif


class Hematite0001Surface(mb.Compound):

    def __init__(self, n_tiles=None, lengths=None, **kwargs):
        super(Hematite0001Surface, self).__init__(**kwargs)

        # load in data from cif file
        cif_filepath = os.path.join(resource_filename("mbuildhelper", "lib"), "surfaces", "fe2o3.cif")
        fe2o3_lattice = load_cif(cif_filepath)

        # check that either n_tiles or lengths are specified
        if n_tiles is None and lengths is None:
            raise ValueError("n_tiles or lengths must be specified")
        if n_tiles is not None and lengths is not None:
            raise ValueError("only one of n_tiles or lengths can be specified")

        # if lengths is specified, determine the number of tiles using the unit cell dimensions
        if lengths is not None:
            _lattice_spacing = np.diagonal(fe2o3_lattice.lattice_vectors) / 10.
            n_tiles = [int(_l / _lattice_spacing[i]) + 1 for i, _l in enumerate(lengths)]

        # number of tiles in y direction must be even so that unit cells line up
        if n_tiles[1] % 2 != 0:
            n_tiles[1] += 1

        # create atoms
        fe = mb.Compound(name='Fe', element='Fe')
        o = mb.Compound(name='O', element='O')

        # populate lattice
        compound_dict = {'O': o, 'Fe': fe}
        fe2o3_lattice = fe2o3_lattice.populate(compound_dict=compound_dict, x=n_tiles[0], y=n_tiles[1], z=n_tiles[2])

        # modify box such that it is rectangular in xy plane
        new_box = mb.Box(lengths=np.diagonal(fe2o3_lattice.box.vectors))
        fe2o3_lattice.box = new_box

        # shift atom positions into new box
        fe2o3_lattice.xyz -= new_box.lengths * np.round_(fe2o3_lattice.xyz / new_box.lengths - 0.5)

        # compute square pair distances
        xyz = fe2o3_lattice.xyz
        rij = xyz - xyz[:, None, :]
        lengths = np.array(new_box.lengths)
        rij = rij - np.array([1, 1, 0]) * lengths * np.round(rij / lengths)
        rij_sqd = np.sum(rij**2, axis=2)

        # determine which pairs are within 3 angstrom cutoff
        rij_sqd[np.tril_indices_from(rij_sqd)] = np.inf # fill lower triangle of square dist matrix with inf to prevent double counting pairs
        pairs = np.transpose(np.where(rij_sqd < .3**2))

        # add bonds for all pairs within cutoff
        for p in pairs:
            fe2o3_lattice.add_bond((fe2o3_lattice[f"Compound[{p[0]}]"], fe2o3_lattice[f"Compound[{p[1]}]"]))

        # add lattice
        self.add(fe2o3_lattice)


if __name__ == '__main__':
    hematite = Hematite0001Surface(lengths=[2, 2, 1])
    hematite.save("hematite_0001.pdb", overwrite=True)
