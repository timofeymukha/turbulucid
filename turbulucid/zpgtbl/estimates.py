from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

__all__ = ["cf_from_rex", "cf_from_redelta99", "redelta99_from_rex",
           "retau_from_rex", "retheta_from_rex", "cf_from_retheta"]


           
# From reX
def cf_from_rex(reX):
    return 0.0283*pow(reX, -2./13)


def redelta99_from_rex(reX):
    return 0.1222*pow(reX, 0.8628)


def retau_from_rex(reX):
    return 0.0145*pow(reX, 0.7858)


def retheta_from_rex(reX):
    return 0.0167*pow(reX, 0.8460)+373.83


# From reDelta99
def cf_from_redelta99(reDelta99):
    return 0.01947*pow(reDelta99, -0.1785)


def rex_from_redelta99(reDelta99):
    return 11.4379*pow(reDelta99, 1.1593)


def retau_from_redelta99(reDelta99):
    return 0.0904*pow(reDelta99, 0.8706)


def retheta_from_redelta99(reDelta99):
    return 0.1312*pow(reDelta99, 0.9808)


# From reTau
def cf_from_retau(reTau):
    return 0.01234*pow(reTau, -0.1960)


def rex_from_retau(reTau):
    return 218.6864*pow(reTau, 1.2726)


def redelta99_from_retau(reTau):
    return 12.7468*pow(reTau, 1.0977)


def retheta_from_retau(reTau):
    return 1.5930*pow(reTau, 1.0766)


# From reTheta
def cf_from_retheta(reTheta):
    return 0.0134*pow(reTheta-373.83, -2./11)


def rex_from_retheta(reTheta):
    return 126.1270*pow(reTheta, 1.1820)


def redelta99_from_retheta(reTheta):
    return 7.9292*pow(reTheta, 1.0196)


def retau_from_retheta(reTheta):
    return 0.5484*pow(reTheta, 0.8877)

