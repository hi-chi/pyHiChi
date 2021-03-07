import sys
sys.path.append("../bin")
import pyHiChi as hichi
import numpy as np
from numba import cfunc, types, carray


min_coords = hichi.Vector3d(-10, -10, 0.0)
max_coords = hichi.Vector3d(10, 10, 0.0)

time_step = 0.05/hichi.c   
grid_size = hichi.Vector3d(128, 1, 1)


# ---------------- create original field ----------

def create_field():
    d = (max_coords.x - min_coords.x) / grid_size.x
    grid_step = hichi.Vector3d(d, d, d)

    @cfunc("float64(float64,float64,float64)")  
    def null_value(x, y, z):
        return 0.0
        
    @cfunc("float64(float64,float64,float64)")
    def field_value(x, y, z):
        return np.sin(0.4*np.pi*x)

    field = hichi.PSATDField(grid_size, min_coords, grid_step, time_step)
    field.set_E(null_value.address, field_value.address, null_value.address)
    field.set_B(null_value.address, null_value.address, field_value.address)
    field.convert_fields_poisson_equation()

    return field

field1 = create_field()


# ------------ apply Python function to field ----------

field2 = create_field()

def func_to_apply_py(x, y, z, field_value):
    field_value.E.y *= np.cos(np.pi/20.0*x)**2
    field_value.B.z *= np.cos(np.pi/20.0*x)**2

field2.apply_function(func_to_apply_py)


# ------------ apply pre-compiled function to field ----------

field3 = create_field()

@cfunc("void(float64,float64,float64,CPointer(float64))")
def func_to_apply_c(x, y, z, field_value_):
    field_value = carray(field_value_, (6,))  # field is array of 6 elements
                                  # (Ex, Ey, Ez, Bx, By, Bz)
    field_value[1] *= np.cos(np.pi/20.0*x)**2  # Ey
    field_value[5] *= np.cos(np.pi/20.0*x)**2  # Bz

field3.apply_function(func_to_apply_c.address)


# --------------- show -----------------

import matplotlib.pyplot as plt

N = 128
x = np.linspace(min_coords.x, max_coords.x, N)

def get_fields(field):
    global x, N
    res = np.zeros(shape=(N,))
    for ix in range(N):
        coord_xy = hichi.Vector3d(x[ix], 0.0, 0.0)
        res[ix] = field.get_E(coord_xy).y
    return res
    
fig = plt.figure()

def plot(field, index, title):
    ax = fig.add_subplot(1,3,index)
    ax.plot(x, get_fields(field))
    ax.set_xlim((min_coords.x, max_coords.x))
    ax.set_xlabel("$x$")
    ax.set_ylabel("$E_y$")
    ax.set_title(title)

plot(field1, 1, "Original field")
plot(field2, 2, "Applied Python function")
plot(field3, 3, "Applied pre-compiled function")

fig.tight_layout()
plt.show()
