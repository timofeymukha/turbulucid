from turbulucid.channel_flow import *
import h5py


# Test if we can open the file, try to grab a dataset.
def test_open_dns_schlatter():
    path = dns_lee_moser()
    dbFile = h5py.File(path, 'r')
    dbFile["180"]["eta"]
