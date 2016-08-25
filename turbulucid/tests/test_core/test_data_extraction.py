from turbulucid import *


def test_interplolate_dataset():
    interpolate_dataset(zpgtbl.dns_schlatter(), 833, "yPlus", "uPlus")
