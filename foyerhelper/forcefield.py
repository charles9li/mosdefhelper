"""ForceField class that extends from foyer's ForceField class."""
import glob
import os
import warnings
from pkg_resources import resource_filename

import foyer


class Forcefield(foyer.Forcefield):
    """Extension of Foyer's ForceField class that provides a number of
    additional automatically-supported force fields.
    """

    def __init__(self, forcefield_files=None, name=None, validation=True, debug=False):
        # first try to use the original foyer Forcefield class
        try:
            super(Forcefield, self).__init__(forcefield_files=forcefield_files,
                                             name=name,
                                             validation=validation,
                                             debug=debug)

        # if original foyer Forcefield class doesn't work, try to search for
        # additional force field files in foyerhelper
        except (FileNotFoundError, OSError):
            all_files_to_load = []
            if forcefield_files is not None:
                if isinstance(forcefield_files, (list, tuple, set)):
                    for file in forcefield_files:
                        all_files_to_load.append(file)
                else:
                    all_files_to_load.append(forcefield_files)

            if name is not None:
                try:
                    file = self.included_forcefields[name]
                except KeyError:
                    raise IOError(f"Forcefield {name} cannot be found")
                else:
                    all_files_to_load.append(file)

            super(Forcefield, self).__init__(forcefield_files=all_files_to_load,
                                             validation=validation,
                                             debug=debug)

    @property
    def included_forcefields(self):
        """Overrides parent property to include additional force fields in
        foyerhelper.
        """
        if len(self._included_forcefields) > 2:
            return self._included_forcefields

        self._included_forcefields = super(Forcefield, self).included_forcefields

        ff_dir = resource_filename("foyerhelper", "forcefields")
        ff_filepaths = set(glob.glob(os.path.join(ff_dir, "xml/*.xml")))

        for ff_filepath in ff_filepaths:
            _, ff_file = os.path.split(ff_filepath)
            basename, _ = os.path.splitext(ff_file)
            self._included_forcefields[basename] = ff_filepath
        return self._included_forcefields

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
            structure = super(Forcefield, self).apply(structure,
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


# different spelling
ForceField = Forcefield
