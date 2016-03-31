import os
import sys

__all__ = ["dns_schlatter"]

def dns_schlatter():

    return os.path.join(sys.modules[__name__].__file__, "..","datasets",
                        "zpgtbl_dns_schlatter.hdf5")
