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
                 edge_smoothing_angle=0.1,
                 wavelength=1.0,
                 total_power=1.0):
        
        self.wavelength = wavelength
        self.pulselength = pulselength*wavelength
        self.phase = phase
        self.R0 = R0*wavelength
        self.total_power = total_power
        self.f_number = f_number
        self.edge_smoothing_angle = edge_smoothing_angle
        self.opening_angle = np.arctan(1.0/(2.0*f_number))
        
        opening_angle = self.opening_angle
        
        @njit
        def longitudinal_field_variation(x):
            return np.sin(2*hichi.pi*x/wavelength + phase) * \
                np.cos(hichi.pi*x/pulselength)**2 * \
                hp.block(x, -0.5*pulselength, 0.5*pulselength)
        
        @njit
        def transverse_shape(angle):
            angle1 = opening_angle - edge_smoothing_angle*0.5
            angle2 = opening_angle + edge_smoothing_angle*0.5
            return hp.block(angle, -1.0, angle1) + \
                (np.cos(0.5*hichi.pi*(angle - angle1)/edge_smoothing_angle)**2 * \
                 hp.block(angle, angle1, angle2) if (not edge_smoothing_angle == 0.0) \
                 else 0.0)
        
        amp = np.sqrt(total_power*4.0/(1.0 - np.cos(opening_angle)))
        
        @njit
        def get_polarisation():
            return hp.Vector3d(0.0, 1.0, 0.0)
               
        @njit
        def mask(x, y, z):
            R = np.sqrt(x*x + y*y + z*z)
            if(R > 1e-5):
                angle = np.arcsin(np.sqrt(y*y + z*z)/R)
                return (amp/R)*longitudinal_field_variation(R - R0)*transverse_shape(angle)*(x < 0)
            else:
                return 0
          
        @cfunc("void(float64,float64,float64,types.CPointer(float64))", nopython=True)
        def EM_field(x, y, z, field_arr_):
            r = hp.Vector3d(x, y, z)
            s1 = hp.cross(get_polarisation(), r)
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
        
        def set_field_(field):
            field.set(EM_field.address)
                
        self.set_field = set_field_
    
    
# --------- creating spherical pulse from C++ ---------------------------------------------------------------


class SphericalPulseC():
    
    def __init__(self, f_number=0.3, R0=16,
                 pulselength=2.0, phase=0.0,
                 edge_smoothing_angle=0.1,
                 wavelength=1.0,
                 total_power=hichi.c):
        
        self.wavelength = wavelength
        self.pulselength = pulselength*wavelength
        self.phase = phase
        self.R0 = R0*wavelength
        self.total_power = total_power
        self.f_number = f_number
        self.edge_smoothing_angle = edge_smoothing_angle
        self.opening_angle = np.arctan(1.0/(2.0*f_number))
        
        tight_focusing_conf = hichi.TightFocusingField(f_number,
                                                       R0,
                                                       wavelength,
                                                       pulselength,
                                                       phase,
                                                       total_power,
                                                       edge_smoothing_angle
                                                      ) 
        def set_field_(field):
            field.set(tight_focusing_conf)
            
        self.set_field = set_field_
