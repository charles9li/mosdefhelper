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
__all__ = ['convert_units', 'convert_torsion']

import numpy as np
from scipy.constants import R
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


# Conversion matrices
# =============================================================================
#   Each function in this section takes in the number of basis functions used
#   and returns a conversion matrix. The functions are added to the appropriate
#   edge in the torsion form graph.


def _opls_to_cosine(n_bases):
    """Returns conversion matrix for OPLS to cosine torsion.

    Parameters
    ----------
    n_bases : int
        Number of basis functions used.

    Returns
    -------
    matrix : np.ndarray, shape (n_bases, n_bases)
        Transformation matrix.
    """
    # initialize transformation matrix
    matrix = np.zeros((n_bases, n_bases))

    # fill in matrix
    for i in range(n_bases):
        matrix[0, i] = 1
        if i > 0:
            matrix[i, i] = (-1)**(i+1)

    return matrix


def _cosine_to_opls(n_bases):
    """Returns conversion matrix for cosine to OPLS torsion.

    Parameters
    ----------
    n_bases : int
        Number of basis functions used.

    Returns
    -------
    matrix : np.ndarray, shape (n_bases, n_bases)
        Transformation matrix.
    """
    return np.linalg.inv(_opls_to_cosine(n_bases))


TORSION_FORM_GRAPH.add_edge('opls', 'cosine', function=_opls_to_cosine)
TORSION_FORM_GRAPH.add_edge('cosine', 'opls', function=_cosine_to_opls)


def _cosine_to_power(n_bases, inv=False):
    """Returns conversion matrix for cosine to power torsion.

    Coefficients of Chebyshev polynomials of the first kind are used to
    construct the transformation matrix due to the polynomials having the
    following property:

        T_n(cos(x)) = cos(n * x)

    This allows for mappings between cosine series into a power series.

    Parameters
    ----------
    n_bases : int
        Number of basis functions used.

    Returns
    -------
    matrix : np.ndarray, shape (n_bases, n_bases)
        Transformation matrix.
    """
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


def _power_to_cosine(n_bases):
    """Returns conversion matrix for cosine to power torsion.

    Parameters
    ----------
    n_bases : int
        Number of basis functions used.

    Returns
    -------
    matrix : np.ndarray, shape (n_bases, n_bases)
        Transformation matrix.
    """
    return np.linalg.inv(_cosine_to_opls(n_bases))


TORSION_FORM_GRAPH.add_edge('cosine', 'power', function=_cosine_to_power)
TORSION_FORM_GRAPH.add_edge('power', 'cosine', function=_power_to_cosine)


def _power_to_rb(n_bases):
    """Returns conversion matrix for power to Ryckaert-Bellemans torsion.

    Note that the reverse conversion uses the exact same matrix.

    Parameters
    ----------
    n_bases : int
        Number of basis functions used.

    Returns
    -------
    matrix : np.ndarray, shape (n_bases, n_bases)
        Transformation matrix.
    """
    matrix = np.identity(n_bases)
    matrix[1::2, 1::2] *= -1
    return matrix


TORSION_FORM_GRAPH.add_edge('power', 'rb', function=_power_to_rb)
TORSION_FORM_GRAPH.add_edge('rb', 'power', function=_power_to_rb)


# Unit conversions
# =============================================================================
#   The unit conversion functions similarly to the torsion conversion, where a
#   units make up the nodes of a graph and edges contain functions that convert
#   between units. The supported units are
#    - kJ/mol
#    - kcal/mol
#    - kB
#   The kB units refer to when reported coefficients are normalized by kB,
#   i.e. cn/kB, which have units of kelvin.

# supported units
UNITS = ['kj/mol', 'kcal/mol', 'kb']
# stylized names used for error checking
UNIT_TO_FULL_NAME = {'kj/mol': 'kJ/mol',
                     'kcal/mol': 'kcal/mol',
                     'kb': 'kB'}

# create graph of units and add nodes and conversion functions
UNITS_GRAPH = nx.DiGraph()
for u in UNITS:
    UNITS_GRAPH.add_node(u)
UNITS_GRAPH.add_edge('kj/mol', 'kcal/mol', function=lambda x: x / 4.184)
UNITS_GRAPH.add_edge('kcal/mol', 'kj/mol', function=lambda x: x * 4.184)
UNITS_GRAPH.add_edge('kb', 'kj/mol', function=lambda x: x * R / 1000)
UNITS_GRAPH.add_edge('kj/mol', 'kb', function=lambda x: x / R * 1000)


def convert_units(x, input_units, output_units):
    """Performs unit conversions on a provided quantity.

    Parameters
    ----------
    x : int or float or array_like
        Quantity that will be scaled for a unit conversion.
    input_units, output_units : str
        Current and desired units of `x`, respectively. One of the following
        string values.

        'kj/mol'
            kJ/mol
        'kcal/mol'
            kcal/mol
        'kb'
            kB units

    Returns
    -------
    x : int or float or array_like
        Quantity with units specified using `output_units`.
    """
    # set to lowercase
    input_units = input_units.lower()
    output_units = output_units.lower()

    # check that input_units and output_units are proper
    if input_units not in UNITS:
        raise ValueError(f"{input_units} is not a valid torsion form")
    if output_units not in UNITS:
        raise ValueError(f"{output_units} is not a valid torsion form")

    # try to find if series of conversions exists between input units and output units
    try:
        conversion_path = nx.shortest_path(UNITS_GRAPH, source=input_units, target=output_units)
    except nx.exception.NetworkXNoPath:
        raise ValueError(f"no conversion method from {UNIT_TO_FULL_NAME[input_units]} torsion "
                         f"to {UNIT_TO_FULL_NAME[output_units]} torsion exists")

    # carry out conversion strategy by traversing conversion path
    for prev_form, curr_form in zip(conversion_path[:-1], conversion_path[1:]):
        x = UNITS_GRAPH[prev_form][curr_form]['function'](x)

    return x


# Main conversion function
# ========================


def convert_torsion(coefficients, input_form, output_form, input_units=None, output_units=None):
    """Takes coefficients for a specified torsion and applies conversion(s) to
    transform them into coefficients for another torsion form.

    Parameters
    ----------
    coefficients : array_like, shape (N,)
        Coefficients corresponding to the torsion form specified by `input_form`.
    input_form, output_form : str
        Torsion form of the input and output coefficients, respectively. One of
        the following string values.

        'opls'
            OPLS torsion form
        'cosine'
            cosine torsion form
        'power'
            power torsion form, consists of linear combination of cosine powers
        'rb'
            Ryckaert-Bellemans torsion form
    input_units, output_units : str, optional
        Current and desired units of `x`, respectively. One of the following
        string values.

        'kj/mol'
            kJ/mol
        'kcal/mol'
            kcal/mol
        'kb'
            kB units

    Returns
    -------
    converted_coefficients : array_like, shape (N,)
        Coefficients corresponding to the torsion form specified by `input_form`.

    Raises
    ------
    ValueError
        If `input_form` or `output_form` aren't valid torsion forms or if a conversion
        strategy can't be found for `input_form` to t`output_form.`
    """
    # set to lowercase
    input_form = input_form.lower()
    output_form = output_form.lower()

    # check that input_form and output_form are proper
    if input_form not in TORSION_FORMS:
        raise ValueError(f"{input_form} is not a valid torsion form")
    if output_form not in TORSION_FORMS:
        raise ValueError(f"{output_form} is not a valid torsion form")

    # determine number of basis functions used
    n_bases = len(coefficients)

    # try to find if series of conversions exists between input_form and output_form
    try:
        conversion_path = nx.shortest_path(TORSION_FORM_GRAPH, source=input_form, target=output_form)
    except nx.exception.NetworkXNoPath:
        raise ValueError(f"no conversion method from {TORSION_FORM_TO_FULL_NAME[input_form]} torsion "
                         f"to {TORSION_FORM_TO_FULL_NAME[output_form]} torsion exists")

    # carry out conversion strategy by traversing conversion path
    for prev_form, curr_form in zip(conversion_path[:-1], conversion_path[1:]):
        matrix = TORSION_FORM_GRAPH[prev_form][curr_form]['function'](n_bases)
        coefficients = np.dot(matrix, coefficients)

    # do unit conversion if specified
    if input_units is not None or output_units is not None:
        coefficients = convert_units(coefficients, input_units, output_units)

    return coefficients


if __name__ == '__main__':
    # print(convert_torsion([0.0, 355.03, -68.19, 791.32], 'opls', 'rb'))
    # print(convert_torsion([-251.06, 428.73, -111.85, 441.27], 'opls', 'rb'))
    # print(convert_torsion([1820.74, -414.41, -1373.14, -30.19, 0.0], 'cosine', 'rb'))
    # print(convert_torsion([2029.99, -751.83, -538.95, -22.10, -51.27], 'cosine', 'rb'))
    # print(convert_torsion([893.21, 176.62, 53.34, 769.93, 0.0], 'cosine', 'rb'))
    # print(convert_torsion([2035.58, -736.9, 57.84, -293.23], 'opls', 'rb'))
    # print(convert_torsion([0.0, 2158.0, 2098.0, 197.3], 'opls', 'rb'))
    print(convert_torsion(np.array([0.0, 710.06, -136.38, 1582.64])/2, 'opls', 'rb', input_units='kb', output_units='kj/mol'))
    print(convert_torsion([630.0+1562.4, -630.0, -1562.4], 'power', 'rb', input_units='kb', output_units='kj/mol'))
    print(convert_torsion([630.0+1562.4, 630.0, -1562.4], 'power', 'rb', input_units='kb', output_units='kj/mol'))
