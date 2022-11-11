"""ForceField class that extends from foyer's ForceField class."""
import os
import warnings
from pkg_resources import resource_filename

import foyer


class ForceField(foyer.Forcefield):
    """Extension of Foyer's ForceField class that provides a number of
    additional automatically-supported force fields.
    """

    def __init__(self, forcefield_files=None, name=None, validation=True, debug=False):
        try:
            super(ForceField, self).__init__(forcefield_files=forcefield_files, name=name, validation=validation, debug=debug)
        except (FileNotFoundError, OSError) as err:
            ff_filepath = resource_filename("foyerhelper", "forcefields")
            if forcefield_files is not None:
                if isinstance(forcefield_files, str):
                    forcefield_files = [forcefield_files]
                forcefield_files = [os.path.join(ff_filepath, "xml", fff) for fff in forcefield_files]
                super(ForceField, self).__init__(forcefield_files=forcefield_files)
            elif name is not None:
                super(ForceField, self).__init__(forcefield_files=os.path.join(ff_filepath, "xml", name))
            else:
                raise err

    def apply(
        self,
        structure,
        references_file=None,
        use_residue_map=True,
        assert_bond_params=True,
        assert_angle_params=True,
        assert_dihedral_params=True,
        assert_improper_params=False,
        verbose=False,
        set_masses_to_zero=False,
        *args,
        **kwargs,
    ):
        """Apply the force field to a molecular structure. Same as
        foyer.ForceField's apply method, but contains an additional keyword
        argument that sets all the masses to zero.
        """
        with warnings.catch_warnings():
            if not verbose:
                warnings.simplefilter("ignore")
            structure = super(ForceField, self).apply(structure,
                                                      references_file=references_file,
                                                      use_residue_map=use_residue_map,
                                                      assert_bond_params=assert_bond_params,
                                                      assert_angle_params=assert_angle_params,
                                                      assert_dihedral_params=assert_dihedral_params,
                                                      assert_improper_params=assert_improper_params,
                                                      verbose=False,
                                                      *args,
                                                      **kwargs)

        # set masses to zero if specified
        if set_masses_to_zero:
            for atom in structure.atoms:
                atom.mass = 0.0

        return structure
