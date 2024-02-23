"""Module for creating SimpleTiledCompounds."""
__all__ = ['SimpleTiledCompound']

import itertools as it

import numpy as np

import mbuild as mb


class SimpleTiledCompound(mb.Compound):
    """Replicates a Compound in any cartesian direction(s).

    This is a simpler version of the TiledCompound recipe in mBuild where
    connectivity between tiles and periodicity is NOT accounted for.

    Parameters
    ----------
    tile : mb.Compound
        The Compound to be replicated.
    n_tiles : array-like, shape=(3,), dtype=int, optional, default=(1, 1, 1)
        Number of times to replicate tile in the x, y and z-directions.
    name : str, optional, default=tile.name
        Descriptive string for the compound.
    """

    def __init__(self, tile, n_tiles, name=None, **kwargs):
        super(SimpleTiledCompound, self).__init__(**kwargs)

        n_tiles = np.asarray(n_tiles)
        periodicity = np.asarray(tile.periodicity)
        if not np.all(n_tiles > 0):
            raise ValueError("Number of tiles must be positive.")

        if tile.box is None:
            raise ValueError(
                "Tile must have a specified box"
            )

        if name is None:
            name = tile.name + "-".join(str(d) for d in n_tiles)
        self.name = name
        self.periodicity = tile.periodicity
        self.box = mb.Box(
            np.array(tile.box.lengths) * n_tiles, angles=tile.box.angles
        )

        # For every tile, assign temporary ID's to particles which are internal
        # to that tile. E.g., when replicating a tile with 1800 particles, every
        # tile will contain particles with ID's from 0-1799. These ID's are used
        # below to fix bonds crossing periodic boundary conditions where a new
        # tile has been placed.
        for idx, particle in enumerate(tile.particles(include_ports=True)):
            particle.index = idx

        # Replicate and place periodic tiles.
        # -----------------------------------
        for ijk in it.product(*[range(i) for i in n_tiles]):
            new_tile = mb.clone(tile)
            new_tile.translate(np.multiply(ijk, np.asarray(tile.box.lengths)))

            self._add_tile(new_tile, ijk)

    def _add_tile(self, new_tile, ijk):
        """Add a tile with a label indicating its tiling position."""
        i, j, k = ijk
        tile_label = f"{self.name}_{i}-{j}-{k}"
        self.add(new_tile, label=tile_label, inherit_periodicity=False)
