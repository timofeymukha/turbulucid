import numpy as np
import pytest

from turbulucid.core.quantities import delta_99, delta_star, momentum_thickness


def test_integral_thicknesses_for_linear_profile():
    y = np.linspace(0.0, 1.0, 101)
    velocity = y.copy()

    assert delta_star(y, velocity, u0=1.0) == pytest.approx(0.5)
    assert momentum_thickness(y, velocity, u0=1.0) == pytest.approx(1 / 6)


def test_momentum_thickness_includes_maximum_velocity_sample():
    y = np.array([0.0, 0.5, 1.0, 1.5])
    velocity = np.array([0.0, 0.5, 1.0, 0.5])

    assert momentum_thickness(y, velocity, u0="max") == pytest.approx(1 / 6)


def test_delta_99_with_interpolation():
    y = np.linspace(0.0, 1.0, 11)
    velocity = y.copy()

    assert delta_99(y, velocity, u0=1.0, interpolate=True) == pytest.approx(
        0.99, abs=2e-4
    )


def test_momentum_thickness_rejects_invalid_cutoff():
    y = np.linspace(0.0, 1.0, 5)

    with pytest.raises(ValueError, match="cutoff"):
        momentum_thickness(y, y, cutoff=5)
