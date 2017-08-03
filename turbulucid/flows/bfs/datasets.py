from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import os
import sys

__all__ = ["jovic_driver", "dns_le_moin"]


def jovic_driver():
    """
    Return the path to the experimental data by Jovic and Driver for the BFS.

    References
    ----------
        Jovic, S., & Driver, D. M. (1994).
        Backward-facing step measurements at low Reynolds number, Re_h= 5000.
        NASA TM 108807,
        1994

        Jovic, S. & Driver, D. (1995).
        Reynolds number effect on the skin friction in separated
        flows behind a backward-Facing Step.
        Experiments in Fluids,
        18(6),
        pp. 464-467


    """

    return os.path.join(os.path.dirname(sys.modules[__name__].__file__),
                        "datasets", "jovic_driver.hdf5")


def dns_le_moin():
    """
    Return the path to the DNS data by Le and Moin for the BFS.

    As found in the ERCOFTAC Classic Collection, case

    References
    ----------
        Le, H. & Moin, P. (1992).
        Direct numerical simulation of turbulent flow over a
        backward-facing step.
        Stanford Univ.,
        Center for Turbulence Research,
        Annual Research Briefs,
        pp. 161-173

    """
    return os.path.join(os.path.dirname(sys.modules[__name__].__file__),
                        "datasets", "dns_le_moin.hdf5")
