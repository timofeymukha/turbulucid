__all__ = ["cf_from_rex", "cf_from_redelta99", "redelta99_from_rex",
           "retau_from_rex", "retheta_from_rex", "cf_from_retheta"]

           
def cf_from_rex(ReX):
    return 0.0282*pow(ReX, -2.0/13)

def cf_from_redelta99(ReDelta99):
    return 0.0203*pow(ReDelta99, -2.0/11)

def redelta99_from_rex(ReX):
    return 0.1635*pow(ReX, 11.0/13)

def retau_from_rex(ReX):
    return 0.0194*pow(ReX, 10.0/13)

def retheta_from_rex(ReX):
    return 0.0167*pow(ReX, 11.0/13)+378.43

def cf_from_retheta(ReTheta):
    return 0.0134*pow(ReTheta - 378.43, -2./11)
