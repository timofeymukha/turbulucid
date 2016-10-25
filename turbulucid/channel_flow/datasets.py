from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os
import sys
import numpy as np

__all__ = ["dns_lee_moser", "dns_lee_moser_retau", "dns_lee_moser_redelta99",
           "dns_lee_moser_redeltastar", "dns_lee_moser_retheta",
           "dns_moser_kim_mansour", "dns_moser_kim_mansour_retau",
           "dns_jimenez", "dns_moser_jimenez_retau"]


def dns_lee_moser():
    """
    Returns the path to the DNS dataset by Lee and Moser for channel flow.
    """

    return os.path.join(sys.modules[__name__].__file__, "..","datasets",
                        "channel_dns_lee_moser.hdf5")


def dns_lee_moser_retau():
    return np.array([182.088, 543.396, 1000.512, 1994.756, 5185.879])


def dns_lee_moser_retheta():
    return np.array([290.308, 1005.680, 1930.604, 3968.745, 10182.024])


def dns_lee_moser_redelta99():
    return np.array([2797.662, 9567.261, 18895.102, 40169.6495, 111714.806])


def dns_lee_moser_redeltastar():
    return np.array([469.305, 1413.505, 2604.424, 5185.416, 12816.696])


def dns_moser_kim_mansour():
    """
    Returns the path to the DNS dataset by Moser, Kim and Mansour for
    channel flow.
    """

    return os.path.join(sys.modules[__name__].__file__, "..","datasets",
                        "channel_moser_kim_mansour.hdf5")


def dns_moser_kim_mansour_retau():
    return np.array([178.12, 392.24, 587.19])


def dns_jimenez():
    """
    Returns the path to the DNS dataset by Jimenez and co-workers for
    channel flow.
    """

    return os.path.join(sys.modules[__name__].__file__, "..","datasets",
                        "channel_jimenez.hdf5")


def dns_moser_jimenez_retau():
    return np.array([180, 550, 934, 2003, 4179])
