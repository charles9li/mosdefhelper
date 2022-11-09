"""Contains utility function for converting between different torsion forms.
Supported forms include:

 - OPLS
     - E(phi) = c0
                + c1 * [ 1 + cos(phi) ]
                + c2 * [ 1 - cos(2 * phi) ]
                + c3 * [ 1 + cos(3 * phi) ]
                + ...
                + cn * [ 1 - (-1)^n * cos(n * phi) ]
                + ...
 - cosine
     - E(phi) = c0
                + c1 * cos(phi)
                + c2 * cos(2 * phi)
                + c3 * cos(3 * phi)
                + ...
                + cn * cos(n * phi)
                + ...
 - power
     - E(phi) = c0
                + c1 * cos(phi)
                + c2 * cos^2(phi)
                + c3 * cos^3(phi)
                + ...
                + cn * cos^n(phi)
                + ...
 - Ryckaert-Bellemans
     - E(psi) = c0
                + c1 * cos(psi)
                + c2 * cos^2(psi)
                + c3 * cos^3(psi)
                + ...
                + cn * cos^n(psi)
                + ...

NOTE: Notice that two different dihedral angle conventions are used, phi and
psi:
 - phi is zero when the first and last particles are on the SAME side of the
   bond formed by the middle two particles (cis conformation)
 - psi is zero when the first and last particles are on OPPOSITE sides (trans
   configuration)
This means that
    phi = psi - pi.
Most of the torsion forms use phi, but due to reasons of convention, the
Ryckaert-Bellemans torsion form uses psi. Note that the power and
Ryckaert-Bellemans torsion are almost the same except for the angle convention
used - to convert between the two, simply multiply each coefficient by (-1)^n.
"""
__all__ = ['convert_torsion']

import numpy as np
import networkx as nx


# currently supported torsion forms
TORSION_FORMS = ['opls', 'cosine', 'power', 'rb']
# used for error message
TORSION_FORM_TO_FULL_NAME = {'opls': 'OPLS',
                             'cosine': 'cosine',
                             'power': 'power',
                             'rb': 'Ryckaert-Bellemans'}


# Create a graph where the nodes are the various torsion forms and the edges
# contain functions that return conversion matrices. If a path exists between
# two forms, then a conversion between those forms is possible. A directed
# graph is used to allow for conversions in both directions and support the
# possibility of a reverse transformation not existing.
TORSION_FORM_GRAPH = nx.DiGraph()
for tf in TORSION_FORMS:
    TORSION_FORM_GRAPH.add_node(tf)


# Create functions that return conversion matrices for a specified
# transformation and add them to the edges of the graph
# ================================================================


def _opls_to_cosine(n_bases, inv=False):
    """Returns conversion matrix for OPLS to cosine torsion."""
    # initialize transformation matrix
    matrix = np.zeros((n_bases, n_bases))

    # fill in matrix
    for i in range(n_bases):
        matrix[0, i] = 1
        if i > 0:
            matrix[i, i] = (-1)**(i+1)

    if not inv:
        return matrix
    else:
        return np.linalg.inv(matrix)


TORSION_FORM_GRAPH.add_edge('opls', 'cosine', function=_opls_to_cosine)
TORSION_FORM_GRAPH.add_edge('cosine', 'opls', function=lambda n_bases: _opls_to_cosine(n_bases, inv=True))


def _cosine_to_power(n_bases, inv=False):
    """Returns conversion matrix for cosine to power torsion."""
    # initialize transformation matrix
    matrix = np.zeros((n_bases, n_bases))

    # fill in matrix
    for i in range(n_bases):
        c = np.polynomial.chebyshev.cheb2poly([0]*i + [1])
        matrix[:len(c), i] = c

    if not inv:
        return matrix
    else:
        return np.linalg.inv(matrix)


TORSION_FORM_GRAPH.add_edge('cosine', 'power', function=_cosine_to_power)
TORSION_FORM_GRAPH.add_edge('power', 'cosine', function=lambda n_bases: _cosine_to_power(n_bases, inv=True))


def _power_to_rb(n_bases):
    """Returns conversion matrix for power to Ryckaert-Bellemans torsion.

    Note that the reverse conversion uses the exact same matrix.
    """
    matrix = np.identity(n_bases)
    matrix[1::2, 1::2] *= -1
    return matrix


TORSION_FORM_GRAPH.add_edge('power', 'rb', function=_power_to_rb)
TORSION_FORM_GRAPH.add_edge('rb', 'power', function=_power_to_rb)


# Main conversion function
# ========================


def convert_torsion(coefficients, from_form, to_form):
    """Takes coefficients for a specified torsion and applies conversion(s) to
    transform them into coefficients for another torsion form.

    Parameters
    ----------
    coefficients : array_like, shape (N,)
        Coefficients corresponding to the torsion form specified by from_form.
    from_form : str
        One of the following string values.

        'opls'
            OPLS torsion form
        'cosine'
            cosine torsion form
        'power'
            power torsion form, consists of linear combination of cosine powers
        'rb'
            Ryckaert-Bellemans torsion form
    to_form : str
        Same string values as from_form.

    Returns
    -------
    converted_coefficients : array_like, shape (N,)
        Coefficients corresponding to the torsion form specified by to_form.

    Raises
    ------
    ValueError
        If from_form or to_form aren't valid torsion forms or if a conversion
        strategy can't be found for from_form to to_form.
    """
    # set to lowercase
    from_form = from_form.lower()
    to_form = to_form.lower()

    # check that from_form and to_form are proper
    if from_form not in TORSION_FORMS:
        raise ValueError(f"{from_form} is not a valid torsion form")
    if to_form not in TORSION_FORMS:
        raise ValueError(f"{to_form} is not a valid torsion form")

    # determine number of basis functions used
    n_bases = len(coefficients)

    # try to find if series of conversions exists between from_form and to_form
    try:
        conversion_path = nx.shortest_path(TORSION_FORM_GRAPH, source=from_form, target=to_form)
    except nx.exception.NetworkXNoPath:
        raise ValueError(f"no conversion method from {TORSION_FORM_TO_FULL_NAME[from_form]} torsion "
                         f"to {TORSION_FORM_TO_FULL_NAME[to_form]} torsion exists")

    # carry out conversion strategy by traversing conversion path
    for prev_form, curr_form in zip(conversion_path[:-1], conversion_path[1:]):
        matrix = TORSION_FORM_GRAPH[prev_form][curr_form]['function'](n_bases)
        coefficients = np.dot(matrix, coefficients)

    return coefficients
