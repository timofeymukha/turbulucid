from turbulucid.zpgtbl import *
import h5py

def test_open_dns_schlatter():
    path = dns_schlatter()
    file = h5py.File(path, 'r')
    eta = file["670"]["eta"]
