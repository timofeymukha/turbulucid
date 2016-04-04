import os
import sys
import h5py

__all__ = ["dns_schlatter", "dns_schlatter_retheta", "dns_schlatter_retau",
           "dns_schlatter_redeltastar", "dns_schlatter_cf",
           "dns_schlatter_shape_factor"]


def dns_schlatter():
    """
    Return the path to the DNS dataset by Schlatter et al for the ZPGTBL.
    """

    return os.path.join(sys.modules[__name__].__file__, "..","datasets",
                        "zpgtbl_dns_schlatter.hdf5")


def dns_schlatter_retheta():
    """
    Return the list of momentum thickness-based Reynolds numbers for the
    profiles in the DNS dataset by Schlatter et al for the ZPGTBL.
    """

    groups = ["670", "1000", "1410", "2000", "2540", "3030", "3270", "3630",
              "4060"]
    dbFile = h5py.File(dns_schlatter())

    reTheta = []
    for i in groups:
        reTheta.append(dbFile[i].attrs["reTheta"])

    return reTheta


def dns_schlatter_redeltastar():
    """
    Return the list of displacement thickness-based Reynolds numbers for the
    profiles in the DNS dataset by Schlatter et al for the ZPGTBL.
    """

    groups = ["670", "1000", "1410", "2000", "2540", "3030", "3270", "3630",
              "4060"]
    dbFile = h5py.File(dns_schlatter())

    reDeltaStar = []
    for i in groups:
        reDeltaStar.append(dbFile[i].attrs["reDeltaStar"])

    return reDeltaStar


def dns_schlatter_retau():
    """
    Return the list of friction velocity-based Reynolds numbers for the
    profiles in the DNS dataset by Schlatter et al for the ZPGTBL.
    """

    groups = ["670", "1000", "1410", "2000", "2540", "3030", "3270", "3630",
              "4060"]
    dbFile = h5py.File(dns_schlatter())

    reTau = []
    for i in groups:
        reTau.append(dbFile[i].attrs["reTau"])

    return reTau


def dns_schlatter_cf():
    """
    Return the list of skin-fiction values for the profiles in the DNS
    dataset by Schlatter et al for the ZPGTBL.
    """

    groups = ["670", "1000", "1410", "2000", "2540", "3030", "3270", "3630",
              "4060"]
    dbFile = h5py.File(dns_schlatter())

    cf = []
    for i in groups:
        cf.append(dbFile[i].attrs["cf"])

    return cf


def dns_schlatter_shape_factor():
    """
    Return the list of shape-factor values for the profiles in the DNS
    dataset by Schlatter et al for the ZPGTBL.
    """

    groups = ["670", "1000", "1410", "2000", "2540", "3030", "3270", "3630",
              "4060"]
    dbFile = h5py.File(dns_schlatter())

    H = []
    for i in groups:
        H.append(dbFile[i].attrs["H"])

    return H
