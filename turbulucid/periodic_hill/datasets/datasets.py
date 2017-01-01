from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os
import sys
import numpy as np

__all__ = ["dns_lee_moser"]


def dns_lee_moser():
    """Return the path to the LES dataset by Frohlich et al 2005."""

    return os.path.join(sys.modules[__name__].__file__, "..","datasets",
                        "periodic_hill_les_frochlich.hdf5")


