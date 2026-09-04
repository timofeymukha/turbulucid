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


def test_delta_99_supports_negative_velocity_orientation():
    y = np.linspace(0.0, 1.0, 11)
    velocity = -y

    assert delta_99(y, velocity, u0=-1.0) == pytest.approx(1.0)


def test_momentum_thickness_rejects_invalid_cutoff():
    y = np.linspace(0.0, 1.0, 5)

    with pytest.raises(ValueError, match="cutoff"):
        momentum_thickness(y, y, cutoff=5)


@pytest.mark.parametrize("function", [momentum_thickness, delta_star, delta_99])
@pytest.mark.parametrize(
    ("y", "velocity", "error", "message"),
    [
        ([0.0], [1.0], ValueError, "at least two"),
        ([0.0, 1.0], [0.0], ValueError, "same number"),
        ([[0.0, 1.0]], [[0.0, 1.0]], ValueError, "one-dimensional"),
        ([0.0, 1.0, 0.5], [0.0, 1.0, 0.5], ValueError, "increasing"),
        ([0.0, np.nan], [0.0, 1.0], ValueError, "finite"),
        (["bad", 1.0], [0.0, 1.0], TypeError, "real numbers"),
    ],
)
def test_profile_validation(function, y, velocity, error, message):
    with pytest.raises(error, match=message):
        function(y, velocity)


@pytest.mark.parametrize("function", [momentum_thickness, delta_star, delta_99])
def test_free_stream_velocity_validation(function):
    y = np.linspace(0.0, 1.0, 5)

    with pytest.raises(ValueError, match="u0 must"):
        function(y, y, u0="unknown")
    with pytest.raises(ValueError, match="nonzero"):
        function(y, y, u0=0)
    with pytest.raises(ValueError, match="finite"):
        function(y, y, u0=np.inf)
    with pytest.raises(TypeError, match="numeric scalar"):
        function(y, y, u0=[1])


def test_momentum_thickness_validates_cutoff_type():
    y = np.linspace(0.0, 1.0, 5)

    with pytest.raises(TypeError, match="integer"):
        momentum_thickness(y, y, cutoff=2.5)
    with pytest.raises(ValueError, match="two profile samples"):
        momentum_thickness(y, y, cutoff=0)


def test_profile_functions_validate_interpolation_flag():
    y = np.linspace(0.0, 1.0, 5)

    with pytest.raises(TypeError, match="boolean"):
        delta_star(y, y, interpolate=1)
