"""mBuildHelper utilities for geometrical operations."""
__all__ = ['scale_rij']

import numpy as np


def scale_rij(particle1, particle2, new_dist):
    """Scales the separation vector between two particles to a new length.

    Here, the position vector of particle2 relative to particle1 is scaled to
    a new distance new_rij. Only the position of particle2 will be shifted.
    Note that any other particles attached to particle2 will not have their
    positions adjusted, so proceed with caution.

    Parameters
    ----------
    particle1, particle2 : mb.Particle
        Two particles whose separation vector will be scaled
    new_dist : float
        Magnitude of new separation vector
    """
    rij = particle2.pos - particle1.pos
    rij_new = rij / np.linalg.norm(rij) * new_dist
    rij_diff = rij_new - rij
    particle2.xyz_with_ports += rij_diff
