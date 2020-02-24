import math as ma
import pyHiChi as hichi
from hichi_primitives import *


# Spherical wave can be created by two ways: from Python and from C++
# C++ is faster
# To create spherical wave from Python, call 'createSphericalPulsePython' (defined below)
# To create spherical wave from C++, call 'createSphericalPulseC' (defined below)


# --------- parameters of spherical wave, they are defined by calling one of functions below ------------------------------


wavelength = None  # wavelength = 1 always
pulselength = None
phase = None
R0 = None

F_number = None
edgeSmoothingAngle = None
openingAngle = None

def getDtCGS(dt):  # return dt in CGS  
    return dt*wavelength/hichi.c
    
setField = None  # function that is setting field for grid, it can be defined by calling any function below
                 # this function is main, it should be called from other python scripts


# --------- creating spherical pulse from python -------------------

Ex = None
Ey = None
Ez = None
Bx = None
By = None
Bz = None

def createSphericalPulsePython(F_number_ = 0.3, R0_ = 16, pulselength_ = 2.0, phase_ = 0, edgeSmoothingAngle_ = 0.1):
    
    from numba import cfunc, float64, jit, njit, jitclass
    
    import numpy as np

    global F_number, R0, wavelength, pulselength, phase, totalPower, edgeSmoothingAngle, openingAngle
    global Ex, Ey, Ez, Bx, By, Bz, setField
    
    wavelength = 1.0
    pulselength = pulselength_*wavelength
    phase = phase_
    R0 = R0_*wavelength
    totalPower = hichi.c
    F_number = F_number_
    edgeSmoothingAngle = edgeSmoothingAngle_
    openingAngle = np.arctan(1.0/(2.0*F_number))

    @njit
    def longitudinalFieldVariation(x):
        return np.sin(2*hichi.pi*x/wavelength + phase) * \
            sqr(np.cos(hichi.pi*x/pulselength)) * \
            block(x, -0.5*pulselength, 0.5*pulselength)
    
    @njit
    def transverseShape(angle):
        return block(angle, -1.0, openingAngle - edgeSmoothingAngle*0.5) + \
            (sqr(np.cos(0.5*hichi.pi*(angle - openingAngle + edgeSmoothingAngle*0.5)/edgeSmoothingAngle)) * \
            block(angle, openingAngle - edgeSmoothingAngle*0.5, openingAngle + edgeSmoothingAngle*0.5) \
            if (not edgeSmoothingAngle == 0) else 0.0)
    
    exclusionRadius = wavelength/10.0
    amp = np.sqrt(totalPower*4.0/(hichi.c*(1.0 - np.cos(openingAngle)))) # need to check!!!
    
    @njit
    def getPolarisation():
        return vector3d(0.0, 1.0, 0.0)
           
    @njit
    def mask(x, y, z, t = 0):
        R = np.sqrt(x*x + y*y + z*z)
        if(R > exclusionRadius):
            angle = np.arcsin(np.sqrt(y*y + z*z)/R)
            return (amp/R)*longitudinalFieldVariation(R - R0)*transverseShape(angle)*(x < 0)
        else:
            return 0
      
    @cfunc("float64(float64,float64,float64)", nopython=True)
    def Ex_(x, y, z):
        r = vector3d(x, y, z)
        s1 = cross(getPolarisation(), r)
        s1.normalize()
        s0 = cross(r, s1)
        s0.normalize()
        return mask(x, y, z)*s0.x
        
    @cfunc("float64(float64,float64,float64)", nopython=True)
    def Ey_(x, y, z):
        r = vector3d(x, y, z)
        s1 = cross(getPolarisation(), r)
        s1.normalize()
        s0 = cross(r, s1)
        s0.normalize()
        return mask(x, y, z)*s0.y
        
    @cfunc("float64(float64,float64,float64)", nopython=True)
    def Ez_(x, y, z):
        r = vector3d(x, y, z)
        s1 = cross(getPolarisation(), r)
        s1.normalize()
        s0 = cross(r, s1)
        s0.normalize()
        return mask(x, y, z)*s0.z
        
    @cfunc("float64(float64,float64,float64)", nopython=True)
    def Bx_(x, y, z):
        r = vector3d(x, y, z)
        s1 = cross(getPolarisation(), r)
        s1.normalize()
        return mask(x, y, z)*s1.x
        
    @cfunc("float64(float64,float64,float64)", nopython=True)
    def By_(x, y, z):
        r = vector3d(x, y, z)
        s1 = cross(getPolarisation(), r)
        s1.normalize()
        return mask(x, y, z)*s1.y
        
    @cfunc("float64(float64,float64,float64)", nopython=True)
    def Bz_(x, y, z):
        r = vector3d(x, y, z)
        s1 = cross(getPolarisation(), r)
        s1.normalize()
        return mask(x, y, z)*s1.z
             
    Ex = Ex_.address   
    Ey = Ey_.address    
    Ez = Ez_.address   
    Bx = Bx_.address    
    By = By_.address
    Bz = Bz_.address
    
    def setField_(grid):
        grid.setE(Ex, Ey, Ez)
        grid.setB(Bx, By, Bz)
        
    setField = setField_
    
    
# --------- creating spherical pulse from C++ ---------------------------------------------------------------


def createSphericalPulseC(F_number_ = 0.3, R0_ = 16, pulselength_ = 2.0, phase_ = 0, edgeSmoothingAngle_ = 0.1):

    global F_number, R0, wavelength, pulselength, phase, totalPower, edgeSmoothingAngle, openingAngle
    global setField
    
    wavelength = 1.0
    pulselength = pulselength_*wavelength
    phase = phase_
    R0 = R0_*wavelength
    totalPower = hichi.c
    F_number = F_number_
    edgeSmoothingAngle = edgeSmoothingAngle_
    openingAngle = np.arctan(1.0/(2.0*F_number))
    
    tightFocusing = hichi.TightFocusingField(F_number, R0, wavelength, pulselength, phase, totalPower, edgeSmoothingAngle)
    
    def setField_(grid):
        grid.set(tightFocusing)
        
    setField = setField_
