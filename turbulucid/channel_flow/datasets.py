import os
import sys

__all__ = ["dns_lee_moser"]

def dns_lee_moser():
    """
    Returns the path to the DNS dataset by Lee and Moser for channel flow.
    """

    return os.path.join(sys.modules[__name__].__file__, "..","datasets",
                        "channel_dns_lee_moser.hdf5")
