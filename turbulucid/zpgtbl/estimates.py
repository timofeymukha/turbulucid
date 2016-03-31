__all__ = ["cf_from_rex", "cf_from_redelta99", "redelta99_from_rex",
           "retau_from_rex", "retheta_from_rex", "cf_from_retheta"]

           
def cf_from_rex(reX):
    return 0.0282*pow(reX, -2.0/13)

def cf_from_redelta99(reDelta99):
    return 0.0203*pow(reDelta99, -2.0/11)

def redelta99_from_rex(reX):
    return 0.1635*pow(reX, 11.0/13)

def retau_from_rex(reX):
    return 0.0194*pow(reX, 10.0/13)

def retheta_from_rex(reX):
    return 0.0167*pow(reX, 11.0/13)+378.43

def cf_from_retheta(reTheta):
    return 0.0134*pow(reTheta - 378.43, -2./11)
