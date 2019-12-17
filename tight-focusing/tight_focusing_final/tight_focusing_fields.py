import sys
sys.path.append("./../../build/src/pyHiChi/Release")
import pyHiChi as pfc
import numpy as np
from numba import cfunc, float64, jit, njit, jitclass
from hichi_primitives import *

# constants and units
TW = 1e+12 * 1e+7  # erg/s
GeV = 1.602e-3  # CGS units

# field initializing
wavelength = 1e-4
pulseLength = 2e-4
phase = 0
R0 = 16*wavelength
totalPower = 50*TW  # assuming the uniform energy stream within the OpeningAngle (see below)
                    # this yields amplitude sqrt(TotalPower*4/(LightVelocity*(cos(OpeningAngle) - 1)))
F_number = 0.3
edgeSmoothingAngle = 0.1  # in radians
openingAngle = np.arctan(1.0/(2.0*F_number))
timeFieldInit = -R0/pfc.c

@njit
def longitudinalFieldVariation(x_ct):
    return np.sin(2*pfc.pi*x_ct/wavelength + phase) * \
        sqr(np.cos(pfc.pi*x_ct/pulseLength)) * \
        block(x_ct, -0.5*pulseLength, 0.5*pulseLength)

@njit
def transverseShape(angle):
    return block(angle, -1.0, openingAngle - edgeSmoothingAngle*0.5) + \
        sqr(np.cos(0.5*pfc.pi*(angle - openingAngle + edgeSmoothingAngle*0.5)/edgeSmoothingAngle)) * \
        block(angle, openingAngle - edgeSmoothingAngle*0.5, openingAngle + edgeSmoothingAngle*0.5)

exclusionRadius = 1e-5
amp = np.sqrt(totalPower*4.0/(pfc.c*(1.0 - np.cos(openingAngle)))) # need to check!!!

@njit  # the bug of numba: cannot lower global constant object of user's class
def getPolarisation():
    return vector3d(0.0, 1.0, 0.0)
    

@njit
def mask(x, y, z, t = 0):
    R = np.sqrt(x*x + y*y + z*z)
    if(R > exclusionRadius):
        angle = np.arcsin(np.sqrt(y*y + z*z)/R)
        return (amp/R)*longitudinalFieldVariation(R + pfc.c*(t + timeFieldInit))*transverseShape(angle)*(x < 0)
    else:
        return 0


@cfunc("float64(float64,float64,float64)", nopython=True)
def Ex(x, y, z):
    r = vector3d(x, y, z)
    s1 = cross(getPolarisation(), r)
    s1.normalize()
    s0 = cross(r, s1)
    s0.normalize()
    return mask(x, y, z)*s0.x
    
@cfunc("float64(float64,float64,float64)", nopython=True)
def Ey(x, y, z):
    r = vector3d(x, y, z)
    s1 = cross(getPolarisation(), r)
    s1.normalize()
    s0 = cross(r, s1)
    s0.normalize()
    return mask(x, y, z)*s0.y
    
@cfunc("float64(float64,float64,float64)", nopython=True)
def Ez(x, y, z):
    r = vector3d(x, y, z)
    s1 = cross(getPolarisation(), r)
    s1.normalize()
    s0 = cross(r, s1)
    s0.normalize()
    return mask(x, y, z)*s0.z
    
@cfunc("float64(float64,float64,float64)", nopython=True)
def Bx(x, y, z):
    r = vector3d(x, y, z)
    s1 = cross(getPolarisation(), r)
    s1.normalize()
    return mask(x, y, z)*s1.x
    
@cfunc("float64(float64,float64,float64)", nopython=True)
def By(x, y, z):
    r = vector3d(x, y, z)
    s1 = cross(getPolarisation(), r)
    s1.normalize()
    return mask(x, y, z)*s1.y
    
@cfunc("float64(float64,float64,float64)", nopython=True)
def Bz(x, y, z):
    r = vector3d(x, y, z)
    s1 = cross(getPolarisation(), r)
    s1.normalize()
    return mask(x, y, z)*s1.z
    
def getFieldFuncs():
    return [[Ex.address, Ey.address, Ez.address], [Bx.address, By.address, Bz.address]]
