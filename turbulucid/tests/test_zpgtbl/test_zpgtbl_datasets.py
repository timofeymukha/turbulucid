from turbulucid.zpgtbl import *
import h5py


# Test if we can open the file, try to grab a dataset.
def test_open_dns_schlatter():
    path = dns_schlatter()
    dbFile = h5py.File(path, 'r')
    dbFile["670"]["eta"]
