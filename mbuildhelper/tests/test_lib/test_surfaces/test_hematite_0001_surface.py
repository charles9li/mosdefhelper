import unittest

import warnings

import numpy as np

import mbuild as mb

from mbuildhelper.lib.molecules import DodecaneUA
from mbuildhelper.lib.recipes import PolyAlkylAcrylateUA
from mbuildhelper.lib.surfaces import Hematite0001Surface
import foyerhelper


class TestHematite0001Surface(unittest.TestCase):

    def test_bonded(self):
        # make 4nm x 4nm x 1nm hematite slab
        hematite = Hematite0001Surface(lengths=[4, 4, 1])
        hematite.name = "HEM"

        # make Lz long
        lengths = list(hematite.box.lengths)
        lengths[2] += 5
        hematite.box = mb.Box(lengths=lengths)

        # apply force field
        iron_oxide = foyerhelper.Forcefield(name="iron-oxide")
        hematite_pmd = iron_oxide.apply(hematite,
                                        assert_angle_params=False,
                                        assert_dihedral_params=False,
                                        residues="HEM")

        # save to files
        hematite_pmd.save("hematite.pdb", overwrite=True)
        hematite_pmd.save("hematite.gro", overwrite=True)
        hematite_pmd.save("hematite.top", overwrite=True)

    def test_dodecane_hematite_interface(self):
        # make 4nm x 4nm x 1nm hematite slab
        hematite_surface = Hematite0001Surface(lengths=[4, 4, 1])

        # make dodecane box
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            lengths = (hematite_surface.box.Lx, hematite_surface.box.Ly, 5)
            dodecane_box = mb.Box(lengths=lengths)
            dodecane_box = mb.fill_box(compound=DodecaneUA(), density=750, box=dodecane_box)
            # dodecane_box = mb.fill_box(compound=DodecaneUA(), n_compounds=2, box=dodecane_box)
            dodecane_box.xyz += np.array([0, 0, hematite_surface.box.Lz])   # shift dodecane above surface

        # apply force fields
        hematite_forcefield = foyerhelper.Forcefield(forcefield_files="iron-oxide.xml")
        hematite_surface = hematite_forcefield.apply(hematite_surface,
                                                     assert_angle_params=False,
                                                     assert_dihedral_params=False,
                                                     set_masses_to_zero=True)
        trappeua = foyerhelper.Forcefield(forcefield_files="trappe-ua-acrylates.xml")
        dodecane_box = trappeua.apply(dodecane_box)

        # add boxes
        interface = dodecane_box + hematite_surface
        interface.box = hematite_surface.box + np.array([0, 0, dodecane_box.box[2], 0, 0, 0])

        interface.save("hematite-dodecane.top", overwrite=True)
        interface.save("hematite-dodecane.gro", overwrite=True)

    def test_10A4_iron_oxide_interface(self):
        # make 5nm x 5nm x 2nm hematite slab
        hematite_surface = Hematite0001Surface(lengths=[5, 5, 2])

        # make box with butyl acrylate
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            lengths = (hematite_surface.box.Lx, hematite_surface.box.Ly, 5)
            butyl_acrylate = mb.Box(lengths=lengths)
            butyl_acrylate = mb.fill_box(PolyAlkylAcrylateUA("10*A4"), n_compounds=1, box=butyl_acrylate)
            butyl_acrylate.xyz += np.array([0, 0, hematite_surface.box.Lz])

        # apply force fields
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            hematite_forcefield = foyerhelper.Forcefield(name="iron-oxide")
            hematite_surface = hematite_forcefield.apply(hematite_surface,
                                                         assert_angle_params=False,
                                                         assert_dihedral_params=False,
                                                         set_masses_to_zero=True)
            trappeua = foyerhelper.Forcefield(name="trappe-ua-acrylates")
            butyl_acrylate = trappeua.apply(butyl_acrylate)

        # add boxes
        interface = butyl_acrylate + hematite_surface
        interface.box = hematite_surface.box + np.array([0, 0, butyl_acrylate.box[2], 0, 0, 0])

        # save to files
        interface.save("hematite-10A4.gro", overwrite=True)
        interface.save("hematite-10A4.top", overwrite=True)


if __name__ == '__main__':
    unittest.main()
