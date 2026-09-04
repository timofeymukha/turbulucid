# This file is part of turbulucid
# (c) 2018 Timofey Mukha
# The code is released under the GNU GPL Version 3 licence.
# See LICENCE.txt and the Legal section in the README for more information

import numpy as np
from scipy.integrate import simpson as simps
from scipy.interpolate import interp1d

__all__ = ["momentum_thickness", "delta_star", "delta_99"]

#: Number of samples the profile is resampled to when interpolate is True.
_INTERPOLATION_SAMPLES = 10000


def _validate_profile(y, v, interpolate):
    """Validate and normalize a one-dimensional velocity profile."""
    try:
        y = np.asarray(y, dtype=float)
        v = np.asarray(v, dtype=float)
    except (TypeError, ValueError) as error:
        raise TypeError("y and v must contain real numbers.") from error

    if y.ndim != 1 or v.ndim != 1:
        raise ValueError("y and v must be one-dimensional arrays.")
    if y.size != v.size:
        raise ValueError("y and v must contain the same number of samples.")
    if y.size < 2:
        raise ValueError("A profile must contain at least two samples.")
    if not np.all(np.isfinite(y)) or not np.all(np.isfinite(v)):
        raise ValueError("Profile coordinates and values must be finite.")
    if np.any(np.diff(y) <= 0):
        raise ValueError("Profile coordinates must be strictly increasing.")
    if not isinstance(interpolate, (bool, np.bool_)):
        raise TypeError("interpolate must be a boolean.")

    return y, v


def _free_stream_velocity(v, u0):
    """Resolve and validate the requested free-stream velocity."""
    if isinstance(u0, str):
        if u0 == "last":
            value = v[-1]
        elif u0 == "max":
            value = np.max(v)
        else:
            raise ValueError("u0 must be 'last', 'max', or a numeric scalar.")
    else:
        if not np.isscalar(u0) or isinstance(u0, (complex, np.complexfloating)):
            raise TypeError("u0 must be 'last', 'max', or a numeric scalar.")
        try:
            value = float(u0)
        except (TypeError, ValueError) as error:
            raise TypeError(
                "u0 must be 'last', 'max', or a numeric scalar."
            ) from error

    if not np.isfinite(value):
        raise ValueError("The free-stream velocity must be finite.")
    if value == 0:
        raise ValueError("The free-stream velocity must be nonzero.")
    return value


def momentum_thickness(y, v, u0="last", cutoff=None, interpolate=False):
    """Compute the momentum thickness.

    Parameters
    ----------
    y : ndarray
        The values of the wall-normal coordinate.
    v : ndarray
        The values of the streamwise velocity.
    u0 : {'last', 'max', value}
        How to compute the free stream velocity. Last will lead to
        using the last value in the v array, max will lead to using
        the maximum value.
    cutoff : int, optional
        Index of the final profile sample included in the integral.
    interpolate : bool
        Whether to add new points to the profile using linear
        interpolation. Useful for coarse profiles.

    Returns
    -------
    float
        The value of the momentum thickness.

    Raises
    ------
    TypeError
        If an option or profile value has an incompatible type.
    ValueError
        If the profile is invalid, unsorted, non-finite, or too short.

    """
    y, v = _validate_profile(y, v, interpolate)
    u0Val = _free_stream_velocity(v, u0)

    if isinstance(u0, str) and u0 == "max":
        cutOff = np.argmax(v)
    else:
        cutOff = v.size - 1

    if cutoff is not None:
        if isinstance(cutoff, (bool, np.bool_)) or not isinstance(
                cutoff, (int, np.integer)):
            raise TypeError("cutoff must be an integer index.")
        cutOff = int(cutoff)

    if cutOff < 1 or cutOff >= v.size:
        raise ValueError("cutoff must include at least two profile samples.")

    y = y[:cutOff + 1]
    v = v[:cutOff + 1]

    if interpolate:
        interp = interp1d(y, v, kind='linear')
        y = np.linspace(y[0], y[-1], _INTERPOLATION_SAMPLES)
        v = interp(y)

    return simps(v/u0Val*(1 - v/u0Val), x=y)


def delta_star(y, v, u0="last", interpolate=False):
    """Compute the displacement thickness.

    Parameters
    ----------
    y : ndarray
        The values of the wall-normal coordinate.
    v : ndarray
        The values of the streamwise velocity.
    u0 : {'last', 'max', value}
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

    Raises
    ------
    TypeError
        If an option or profile value has an incompatible type.
    ValueError
        If the profile is invalid, unsorted, non-finite, or too short.

    """
    y, v = _validate_profile(y, v, interpolate)
    u0Val = _free_stream_velocity(v, u0)

    if interpolate:
        interp = interp1d(y, v, kind='linear')
        y = np.linspace(y[0], y[-1], _INTERPOLATION_SAMPLES)
        v = interp(y)

    return simps(1 - v/u0Val, x=y)


def delta_99(y, v, u0="last", interpolate=False):
    """Compute delta_99.

    Parameters
    ----------
    y : ndarray
        The values of the wall-normal coordinate.
    v : ndarray
        The values of the streamwise velocity.
    u0 : {'last', 'max', value}
        How to compute the free stream velocity. Last will lead to
        using the last value in the v array, max will lead to using
        the maximum value.
    interpolate : bool
        Whether to add new points to the profile using linear
        interpolation. Useful for coarse profiles.

    Returns
    -------
    float
        The value of delta_99. This is always one of the y values of the
        profile, i.e. the first one at which the velocity reaches 99% of
        the free stream value; no interpolation between samples is done.
        The accuracy is therefore limited by the resolution of the
        profile, which is what the interpolate option is for.

    Raises
    ------
    TypeError
        If an option or profile value has an incompatible type.
    ValueError
        If the profile is invalid or the computed value is not positive.

    """
    y, v = _validate_profile(y, v, interpolate)
    u0Val = _free_stream_velocity(v, u0)

    if interpolate:
        interp = interp1d(y, v, kind='linear')
        y = np.linspace(y[0], y[-1], _INTERPOLATION_SAMPLES)
        v = interp(y)

    candidates = np.flatnonzero(v/u0Val >= 0.99)
    if candidates.size == 0:
        raise ValueError("The profile does not reach 99% of u0.")

    delta99 = y[candidates[0]]
    if delta99 <= 0:
        raise ValueError("delta_99 is not positive.")

    return delta99
