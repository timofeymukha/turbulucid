from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

__all__ = ["dean_reb_from_retau", "dean_retau_from_reb", "dean_rec_from_reb",
           "dean_reb_from_rec", "dean_rec_from_retau", "dean_retau_from_rec"]


def uc_over_from_retheta(reTheta):
    """
    Get the ratio of the centerline velocity and the bulk velocity using
    a fit to the DNS data of Lee and Moser.
    """
    return 1 + 0.3427*reTheta**(-0.1287)


def delta_over_theta_from_retheta(reTheta):
    """
    Get the ratio of the half-height and the momentum thickness using a
    fit to the DNS data of Lee and Moser.
    """
    return 11. + 2.603e-4*reTheta**(-0.9834)


def dean_retau_from_reb(reB):
    return 0.175*reB**0.875


def dean_rec_from_reb(reB):
    return 1.27*reB**0.988


def dean_reb_from_retau(reTau):
    return 1/0.175*reTau**(1/0.875)


def dean_reb_from_rec(reC):
    return 1/1.27*reC**(1/0.988)


def dean_retau_from_rec(reC):
    reB = dean_reb_from_rec(reC)
    return dean_retau_from_reb(reB)


def dean_rec_from_retau(reTau):
    reB = dean_reb_from_retau(reTau)
    return dean_rec_from_reb(reB)

