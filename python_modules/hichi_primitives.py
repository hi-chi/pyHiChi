import numpy as np
import os
from numba import cfunc, float64, jit, njit
from  numba.experimental import jitclass

@jitclass([('x', float64), ('y', float64), ('z', float64)])
class Vector3d(object):

    def __init__(self, x = 0.0, y = 0.0, z = 0.0):
        self.x = x
        self.y = y
        self.z = z
        
    def norm(self):
    
        return np.sqrt(self.x**2+self.y**2+self.z**2)
    def norm2(self):
        return self.x**2+self.y**2+self.z**2
        
    def normalize(self):
        norm = self.norm()
        if (norm != 0):
            self.x /= norm
            self.y /= norm
            self.z /= norm

@njit
def block(x, xmin, xmax):
    return (np.sign(x - xmin) + np.sign(xmax - x))*0.5

@njit
def segment(x, x1, y1, x2, y2):
    return block(x, x1, x2)*(y1 + (y2 - y1)*(x - x1)/(x2 - x1))

@njit
def minimum(a, b):
    return 0.5*(b*(np.sign(a - b) + 1.0) + a*(np.sign(b - a) + 1.0))

@njit
def cross(v1, v2):
    return Vector3d(v1.y*v2.z - v1.z*v2.y, \
                    v1.z*v2.x - v1.x*v2.z, \
                    v1.x*v2.y - v1.y*v2.x)

@njit
def sqr(x):
    return x*x
      
def create_dir(dir):
    if (os.path.exists(dir)): 
        for (dirpath, dirnames, filenames) in os.walk(dir):
            for file in filenames:
                os.remove(dir + "/" + file)
    else: os.mkdir(dir)

# useful enums
from enum import Enum

class Axis(Enum):
    X = "x"
    Y = "y"
    Z = "z"

class Field(Enum):
    E = "E"
    B = "B"
    J = "J"    

class Plane(Enum):
    XOY = (Axis.X, Axis.Y)
    XOZ = (Axis.X, Axis.Z)
    YOZ = (Axis.Y, Axis.Z)
    
def get_coord_value(vector, axis):
    if axis == Axis.X:
        return vector.x
    elif axis == Axis.Y:
        return vector.y
    elif axis == Axis.Z:
        return vector.z
    else:
        raise ValueException("ERROR: wrong arg in Visual.get_coord_value") 
