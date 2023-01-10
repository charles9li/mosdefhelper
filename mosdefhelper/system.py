__all__ = ['System']

import numpy as np

import mbuild as mb
import foyer


class _Compound(object):
    """Private container object that stores information of compounds in a
    mosdefhelper.System object.

    Attributes
    ----------
    compound : mb.Compound
        mBuild Compound object.
    n : int
        Number of copies of this compound.
    forcefield : foyer.Forcefield
        Force field to apply to compound.
    apply_kwargs : dict
        Dictionary that contains keyword args when applying force field.
    """

    def __init__(self, compound, n, forcefield, apply_kwargs):
        self.compound = compound
        self.n = n
        self.forcefield = forcefield
        self.apply_kwargs = apply_kwargs


class System(object):
    """Container object that stores information about a collection of compounds
    and their overall composition. This class also contains methods for writing
    out to a number of different file formats.

    Attributes
    ----------
    box : mb.Box
        mbuild Box object that stores information about the size and shape of
        the box.
    compounds : list of mb.Compound
        List of compounds in the system.
    compound_names : list of str
        List of compound names in the system.
    n_compounds : list of int
        List of numbers of each compound.
    """

    def __init__(self, box=None):
        self.box = box
        self._compounds = []

    def add_compound(self, compound, n=1, forcefield=None, apply_kwargs=None):
        """Add a mbuild Compound to the System.

        Parameters
        ----------
        compound : mb.Compound
            mBuild Compound object.
        n : int, optional, default=1
            Number of this compound in the system.
        forcefield : foyer.Forcefield, optional
            Force field to apply to the compound.
        apply_kwargs : dict, optional
            Dictionary that contains keyword args when applying force field.
        """
        if forcefield is None:
            forcefield = foyer.Forcefield()
        if apply_kwargs is None:
            apply_kwargs = dict()
        self._compounds.append(_Compound(compound, n, forcefield, apply_kwargs))

    @property
    def compounds(self):
        return [_c.compound for _c in self._compounds]

    @property
    def compound_names(self):
        return [_c.compound.name for _c in self._compounds]

    @property
    def n_compounds(self):
        return [_c.n for _c in self._compounds]

    def to_parmed(self, shift_positions=None):
        """Creates a ParmEd structure from the information in the system."""
        # fill box
        system = mb.fill_box(compound=self.compounds, n_compounds=self.n_compounds, box=self.box)

        # shift positions if specified
        if shift_positions is not None:
            system.xyz += np.array(shift_positions)

        # separate molecules into different compounds and apply force fields
        child_iter = iter(system.children)
        parametrized_compounds = []
        for _c in self._compounds:
            compound_all_copies = mb.Compound()
            for _ in range(_c.n):
                compound_all_copies.add(mb.clone(next(child_iter)))
            compound_pmd = _c.forcefield.apply(compound_all_copies, **_c.apply_kwargs)
            parametrized_compounds.append(compound_pmd)

        # combine structures
        structure = None
        for pc in parametrized_compounds:
            if structure is None:
                structure = pc
            else:
                structure += pc
        structure.box = np.append(np.array(self.box.lengths) * 10.0, self.box.angles)

        return structure
