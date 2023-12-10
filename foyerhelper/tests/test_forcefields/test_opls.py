import unittest

import mbuild as mb

from mbuildhelper.lib.molecules import StearicAcid
from foyerhelper import Forcefield


class TestOPLS(unittest.TestCase):

    def test_stearic_acid(self):
        stearic_acid = StearicAcid()
        opls = Forcefield(name='oplsaa')
        stearic_acid = opls.apply(stearic_acid)

        print(len(stearic_acid.bond_types))
        bond_types = list(stearic_acid.bond_types)
        for i, bond_type in enumerate(stearic_acid.bond_types):
            print(i, bond_type)
        for bond in stearic_acid.bonds:
            print(bond_types.index(bond.type), bond.atom1, bond.atom2)

        print(len(stearic_acid.angle_types))
        angle_types = list(stearic_acid.angle_types)
        for i, angle_type in enumerate(stearic_acid.angle_types):
            print(i, angle_type)
        for angle in stearic_acid.angles:
            print(angle_types.index(angle.type), angle.atom1, angle.atom2, angle.atom3)

        print(len(stearic_acid.rb_torsion_types))
        rb_torsion_types = list(stearic_acid.rb_torsion_types)
        for i, rb_torsion_type in enumerate(stearic_acid.rb_torsion_types):
            print(i, rb_torsion_type)
        for rb_torsion in stearic_acid.rb_torsions:
            print(rb_torsion_types.index(rb_torsion.type), rb_torsion.atom1, rb_torsion.atom2, rb_torsion.atom3, rb_torsion.atom4)

        for atom in stearic_acid.atoms:
            print(atom.name, atom.sigma, atom.epsilon, atom.charge)

    def test_compare_oplsaa_lopls(self):
        stearic_acid = StearicAcid()
        opls = Forcefield(name="oplsaa")
        lopls = Forcefield(name="l-opls")
        stearic_acid_opls = opls.apply(stearic_acid)
        stearic_acid_lopls = lopls.apply(stearic_acid)

        # make sure bond interactions are equal
        for bond_opls, bond_lopls in zip(stearic_acid_opls.bonds, stearic_acid_lopls.bonds):
            self.assertEqual(bond_opls.type, bond_lopls.type)

        # make sure angle interactions are equal
        for angle_opls, angle_lopls in zip(stearic_acid_opls.angles, stearic_acid_lopls.angles):
            self.assertEqual(angle_opls.type, angle_lopls.type)

        # # make sure torsions are equal
        # for rb_torsion_opls, rb_torsion_lopls in zip(stearic_acid_opls.rb_torsions, stearic_acid_lopls.rb_torsions):
        #     self.assertEqual(rb_torsion_opls.type, rb_torsion_lopls.type)

        for atom_opls, atom_lopls in zip(stearic_acid_opls.atoms, stearic_acid_lopls.atoms):
            print(atom_opls)
            self.assertEqual(atom_opls.sigma, atom_lopls.sigma)
            print(atom_opls.epsilon, atom_lopls.epsilon)
            print(atom_opls.charge, atom_lopls.charge)


if __name__ == '__main__':
    unittest.main()
