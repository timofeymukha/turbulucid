from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from scipy.integrate import simps

__all__ = ["momentum_thickness", "delta_star", "delta_99"]


def momentum_thickness(y, v):
    """Compute the momentum thickness.

    Parameters
    ----------
    y : ndarray
        The values of the wall-normal coordinate.
    v : ndarray
        The values of the streamwise velocity.

    Returns
    -------
    float
        The value of the momentum thickness.

    """
    u0 = v[-1]
    return simps(v/u0*(1-v/u0), x=y)


def delta_star(y, v):
    """Compute the displacement thickness.

    Parameters
    ----------
    y : ndarray
        The values of the wall-normal coordinate.
    v : ndarray
        The values of the streamwise velocity.

    Returns
    -------
    float
        The value of the displacement thickness.

    """
    u0 = v[-1]
    return simps(1-v/u0, x=y)


def delta_99(y, v):
    """Compute delta_99.

    Parameters
    ----------
    y : ndarray
        The values of the wall-normal coordinate.
    v : ndarray
        The values of the streamwise velocity.

    Returns
    -------
    float
        The value of delta_99.

    Raises
    ------
    ValueError
        If the computed value is not positive.

    """
    #interp = interp1d(y, v])
    #newY = np.linspace(y[0], y[-1], 10000)
    #newV = interp(newY)
    u0 = v[-1]
    delta99 = 0
    for i in range(v.size):
        if v[i] >= 0.99*u0:
            delta99 = y[i]
            break

    if delta99 <= 0:
        raise ValueError("delta_99 is not positive!")

    return delta99

