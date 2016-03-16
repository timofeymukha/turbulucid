import numpy as np
from scipy.integrate import simps

__all__ = ["momentum_thickness", "delta_star", "delta_99"]


def momentum_thickness(field, u0, yLoc):
    X = field["X"]
    Y = field["Y"]
    V = field["V"]

    x = X[0,:]
    theta = np.zeros(x.size) 
    for i in xrange(x.size):
        yIdxTop = np.argmin(np.abs(Y[:, i]-yLoc[i]))
        y = Y[0:yIdxTop,i]
        theta[i] = simps(V[0:yIdxTop, i]/u0[i]*(1-V[0:yIdxTop, i]/u0[i]), x=y)     

    return theta   

def delta_star(field, u0, yLoc):
    X = field["X"]
    Y = field["Y"]
    V = field["V"]
    x = X[0,:]

    deltaStar = np.zeros(x.size) 

    for i in xrange(x.size):
        yIdxTop = np.argmin(np.abs(Y[:, i]-yLoc[i]))
        y =Y[0:yIdxTop, i]
        deltaStar[i] = simps((1-V[0:yIdxTop, i]/u0[i]), x=y)

    return deltaStar 

def delta_99(field, u0, yLoc):
    X = field["X"]
    Y = field["Y"]
    V = field["V"]
    x = X[0,:]

    delta99 = np.zeros(x.shape)

    for i in xrange(x.size):
        #interp = interp1d(Y[:,i], V[:,i])
        #newY = np.linspace(Y[0,i], yLoc[i], 10000)
        #newV = interp(newY)
        for j in range(V[:,i].size):
            if V[j, i] >= 0.99*u0[i]:
                delta99[i] = Y[j, i]
                break
    return delta99


