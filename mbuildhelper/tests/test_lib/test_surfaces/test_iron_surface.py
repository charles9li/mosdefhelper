import unittest

import warnings

import numpy as np

import mbuild as mb
import foyer

from mbuildhelper.lib.molecules import DodecaneUA
from mbuildhelper.lib.surfaces import IronSurface
import foyerhelper


class TestIronSurface(unittest.TestCase):

    def test_iron_dodecane_interface(self):
        # make 5nm x 5nm x 1nm iron slab
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            iron_surface = IronSurface(crystal_structure='fcc', temperature=313.15, lengths=[5, 5, 1])

        # make dodecane box
        print("filling dodecane box")
        dodecane_box = mb.Box(lengths=(iron_surface.box.Lx, iron_surface.box.Ly, 5))
        dodecane_box = mb.fill_box(compound=DodecaneUA(), density=750, box=dodecane_box)
        dodecane_box.xyz += np.array([0, 0, iron_surface.box.Lz])   # shift dodecane positions above the iron surface
        print("dodecane_box_filled")

        # apply force fields
        iron_forcefield = foyerhelper.Forcefield(name="fcc-metals")
        iron_surface = iron_forcefield.apply(iron_surface)
        trappeua = foyer.Forcefield(name='trappe-ua')
        dodecane_box = trappeua.apply(dodecane_box)

        # add boxes
        system = iron_surface + dodecane_box
        system.box = iron_surface.box + np.array([0, 0, dodecane_box.box[2], 0, 0, 0])  # replace box with combined box

        system.save('iron-dodecane.pdb', overwrite=True)
        system.save('iron-dodecane.gro', overwrite=True)
        system.save('iron-dodecane.top', overwrite=True)


if __name__ == '__main__':
    unittest.main()
