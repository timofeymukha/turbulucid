# This file is part of turbulucid
# (c) 2018 Timofey Mukha
# The code is released under the GNU GPL Version 3 licence.
# See LICENCE.txt and the Legal section in the README for more information

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from scipy.integrate import simps
from scipy.interpolate import interp1d
import numpy as np

__all__ = ["momentum_thickness", "delta_star", "delta_99"]


def momentum_thickness(y, v, u0="last", interpolate=False):
    """Compute the momentum thickness.

    Parameters
    ----------
    y : ndarray
        The values of the wall-normal coordinate.
    v : ndarray
        The values of the streamwise velocity.
    u0 : {'last', 'max'}
        How to compute the free stream velocity. Last will lead to
        using the last value in the v array, max will lead to using
        the maximum value.
    interpolate : bool
        Whether to add new points to the profile using linear
        interpolation. Useful for coarse profiles.

    Returns
    -------
    float
        The value of the momentum thickness.

    """
    if u0 is "last":
        u0Val = v[-1]
        cutOff = -1
    elif u0 is "max":
        u0Val = np.max(v)
        cutOff = np.argmax(v)
    else:
        raise ValueError("u0 should be either 'last' or 'max'.")

    if interpolate:
        interp = interp1d(y, v, kind='linear')
        y = np.linspace(y[0], y[cutOff], 10000)
        v = interp(y)
    else:
        y = y[:cutOff]
        v = v[:cutOff]

    return simps(v/u0Val*(1-v/u0Val), x=y)


def delta_star(y, v, u0="last", interpolate=False):
    """Compute the displacement thickness.

    Parameters
    ----------
    y : ndarray
        The values of the wall-normal coordinate.
    v : ndarray
        The values of the streamwise velocity.
    u0 : {'last', 'max'}
        How to compute the free stream velocity. Last will lead to
        using the last value in the v array, max will lead to using
        the maximum value.
    interpolate : bool
        Whether to add new points to the profile using linear
        interpolation. Useful for coarse profiles.

    Returns
    -------
    float
        The value of the displacement thickness.

    """
    if u0 is "last":
        u0Val = v[-1]
    elif u0 is "max":
        u0Val = np.max(v)
    else:
        raise ValueError("u0 should be either 'last' or 'max'.")

    if interpolate:
        interp = interp1d(y, v, kind='linear')
        y = np.linspace(y[0], y[-1], 10000)
        v = interp(y)

    return simps(1-v/u0Val, x=y)


def delta_99(y, v, u0="last", interpolate=False):
    """Compute delta_99.

    Parameters
    ----------
    y : ndarray
        The values of the wall-normal coordinate.
    v : ndarray
        The values of the streamwise velocity.
    u0 : {'last', 'max'}
        How to compute the free stream velocity. Last will lead to
        using the last value in the v array, max will lead to using
        the maximum value.
    interpolate : bool
        Whether to add new points to the profile using linear
        interpolation. Useful for coarse profiles.

    Returns
    -------
    float
        The value of delta_99.

    Raises
    ------
    ValueError
        If the computed value is not positive.

    """
    if u0 is "last":
        u0Val = v[-1]
    elif u0 is "max":
        u0Val = np.max(v)
    else:
        raise ValueError("u0 should be either 'last' or 'max'.")

    if interpolate:
        interp = interp1d(y, v, kind='linear')
        y = np.linspace(y[0], y[-1], 10000)
        v = interp(y)

    delta99 = 0
    for i in range(v.size):
        if v[i] >= 0.99*u0Val:
            delta99 = y[i]
            break

    if delta99 <= 0:
        raise ValueError("delta_99 is not positive!")

    return delta99

