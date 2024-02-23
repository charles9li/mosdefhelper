"""Rectilinear version of the Hematite (0001) unit cell."""
__all__ = ['Hematite0001RectilinearUnitCell']

import os
from pkg_resources import resource_filename

import numpy as np

import mbuild as mb
from mbuild.lattice import load_cif


class Hematite0001RectilinearUnitCell(mb.Compound):
    """Rectilinear version of the Hematite (0001) unit cell.

    Contains twice as many atoms as the normal triclinic unit cell.
    """

    def __init__(self, add_bonds=True, **kwargs):
        super(Hematite0001RectilinearUnitCell, self).__init__(**kwargs)

        # load in data from cif file
        cif_filepath = os.path.join(resource_filename("mbuildhelper", "lib"), "surfaces", "fe2o3.cif")
        fe2o3_lattice = load_cif(cif_filepath)

        # create atoms
        fe = mb.Compound(name='Fe', element='Fe')
        o = mb.Compound(name='O', element='O')

        # populate lattice and create unit cell that is replicated once in the y direction
        compound_dict = {'O': o, 'Fe': fe}
        fe2o3_unit_cell = fe2o3_lattice.populate(compound_dict=compound_dict, x=1, y=2, z=1)

        # modify unit cell such that it is rectangular in xy plane
        new_box = mb.Box(lengths=np.diagonal(fe2o3_unit_cell.box.vectors))
        fe2o3_unit_cell.box = new_box

        # shift atom positions into new box
        fe2o3_unit_cell.xyz -= new_box.lengths * np.round_(fe2o3_unit_cell.xyz / new_box.lengths - 0.5)
        fe2o3_unit_cell.periodicity = (True, True, True)

        # add bonds (needed for foyer parameterization)
        if add_bonds:
            prev_particle = None
            for particle in fe2o3_unit_cell.particles():
                if prev_particle is not None:
                    fe2o3_unit_cell.add_bond((prev_particle, particle))
                prev_particle = particle

        # add lattice
        self.add(fe2o3_unit_cell)

        # change name
        self.name = "HEM"


if __name__ == '__main__':
    hematite_unit_cell = Hematite0001RectilinearUnitCell(add_bonds=True)
    hematite_unit_cell.save("hematite_0001_unit_cell.pdb", overwrite=True)
