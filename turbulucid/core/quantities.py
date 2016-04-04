from scipy.integrate import simps

__all__ = ["momentum_thickness", "delta_star", "delta_99"]


def momentum_thickness(y, v):
    u0 = v[-1]
    return simps(v/u0*(1-v/u0), x=y)


def delta_star(y, v):
    u0 = v[-1]
    return simps(1-v/u0, x=y)


def delta_99(y, v):
    #interp = interp1d(y, v])
    #newY = np.linspace(y[0], y[-1], 10000)
    #newV = interp(newY)
    u0 = v[-1]
    delta99 = 0
    for i in range(v.size):
        if v[i] >= 0.99*u0:
            delta99 = y[i]
            break

    assert delta99 > 0

    return delta99

