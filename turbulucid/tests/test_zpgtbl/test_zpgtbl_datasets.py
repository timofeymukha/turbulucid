from turbulucid.zpgtbl import *
import h5py
import numpy as np

# Test if we can open the file, try to grab a dataset.
def test_open_dns_schlatter():
    path = dns_schlatter()
    dbFile = h5py.File(path, 'r')
    dbFile["670"]["eta"]


def test_dns_shlatter_retheta():
    reThetaTrue = np.array([677.452, 1006.534, 1420.960, 2000.703, 2536.827, 3031.843,
                   3273.634, 3626.304, 4061.378])

    reTheta = dns_schlatter_retheta()

    assert np.all(reThetaTrue == reTheta)


def test_dns_shlatter_redeltastar():
    reDeltaStarTrue = np.array([998.106, 1459.397, 2030.877, 2827.938, 3563.325,
                       4237.594, 4567.562, 5044.407, 5633.318])

    reDeltaStar = dns_schlatter_redeltastar()

    assert np.all(reDeltaStarTrue == reDeltaStar)


def test_dns_shlatter_retau():
    reTauTrue = np.array([252.255, 359.3794, 492.2115, 671.124, 830.0115, 974.1849,
                 1043.4272, 1145.1699, 1271.535])

    reTau = dns_schlatter_retau()

    assert np.all(reTauTrue == reTau)

def test_dns_shlatter_cf():
    reCfTrue = np.array([0.004777047, 0.004264437, 0.003884512, 0.003539148,
                0.003327997, 0.003189402, 0.003123259, 0.003055599,
                0.002970989])

    reCf = dns_schlatter_cf()

    assert np.all(reCfTrue == reCf)


def test_dns_shlatter_H():
    reHTrue = np.array([1.473324, 1.449924, 1.429229, 1.413472, 1.404639, 1.397696,
               1.395257, 1.39106, 1.387046])

    reH = dns_schlatter_shape_factor()

    assert np.all(reHTrue == reH)
