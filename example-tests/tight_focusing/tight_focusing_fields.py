import math as ma
import pyHiChi as hichi
import hichi_primitives as hp
from numba import cfunc, float64, jit, njit, jitclass, types, carray 
import numpy as np

# Spherical wave can be created by two ways: with Python and with C++
# C++ is faster

# --------- creating spherical pulse with python -------------------


class SphericalPulsePython():
    
    def __init__(self, f_number=0.3, R0=16,
                 pulselength=2.0, phase=0.0,
                 edgeSmoothingAngle=0.1,
                 wavelength=1.0,
                 totalPower=hichi.c):
        
        self.wavelength = wavelength
        self.pulselength = pulselength*wavelength
        self.phase = phase
        self.R0 = R0*wavelength
        self.totalPower = totalPower
        self.f_number = f_number
        self.edgeSmoothingAngle = edgeSmoothingAngle
        self.openingAngle = np.arctan(1.0/(2.0*f_number))
        
        openingAngle = self.openingAngle
        
        @njit
        def longitudinalFieldVariation(x):
            return np.sin(2*hichi.pi*x/wavelength + phase) * \
                np.cos(hichi.pi*x/pulselength)**2 * \
                hp.block(x, -0.5*pulselength, 0.5*pulselength)
        
        @njit
        def transverseShape(angle):
            angle1 = openingAngle - edgeSmoothingAngle*0.5
            angle2 = openingAngle + edgeSmoothingAngle*0.5
            return hp.block(angle, -1.0, angle1) + \
                (np.cos(0.5*hichi.pi*(angle - angle1)/edgeSmoothingAngle)**2 * \
                 hp.block(angle, angle1, angle2) if (not edgeSmoothingAngle == 0.0) \
                 else 0.0)
        
        amp = np.sqrt(totalPower*4.0/(hichi.c*(1.0 - np.cos(openingAngle))))
        
        @njit
        def getPolarisation():
            return hp.vector3d(0.0, 1.0, 0.0)
               
        @njit
        def mask(x, y, z):
            R = np.sqrt(x*x + y*y + z*z)
            if(R > 1e-5):
                angle = np.arcsin(np.sqrt(y*y + z*z)/R)
                return (amp/R)*longitudinalFieldVariation(R - R0)*transverseShape(angle)*(x < 0)
            else:
                return 0
          
        @cfunc("void(float64,float64,float64,types.CPointer(float64))", nopython=True)
        def EMField(x, y, z, field_arr_):
            r = hp.vector3d(x, y, z)
            s1 = hp.cross(getPolarisation(), r)
            s1.normalize()
            s0 = hp.cross(r, s1)
            s0.normalize()
            m = mask(x, y, z)
            field_arr = carray(field_arr_, (6))         
            field_arr[0] = m*s0.x  # Ex
            field_arr[1] = m*s0.y  # Ey
            field_arr[2] = m*s0.z  # Ez              
            field_arr[3] = m*s1.x  # Bx
            field_arr[4] = m*s1.y  # By
            field_arr[5] = m*s1.z  # Bz
        
        def setField_(grid):
            grid.set(EMField.address)
                
        self.setField = setField_
    
    
# --------- creating spherical pulse from C++ ---------------------------------------------------------------


class SphericalPulseC():
    
    def __init__(self, f_number=0.3, R0=16,
                 pulselength=2.0, phase=0.0,
                 edgeSmoothingAngle=0.1,
                 wavelength=1.0,
                 totalPower=hichi.c):
        
        self.wavelength = wavelength
        self.pulselength = pulselength*wavelength
        self.phase = phase
        self.R0 = R0*wavelength
        self.totalPower = totalPower
        self.f_number = f_number
        self.edgeSmoothingAngle = edgeSmoothingAngle
        self.openingAngle = np.arctan(1.0/(2.0*f_number))
        
        tightFocusing = hichi.TightFocusingField(f_number,
                                                 R0,
                                                 wavelength,
                                                 pulselength,
                                                 phase,
                                                 totalPower,
                                                 edgeSmoothingAngle) 
        def setField_(grid):
            grid.set(tightFocusing)
            
        self.setField = setField_
