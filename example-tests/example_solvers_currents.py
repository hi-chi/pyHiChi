import sys
sys.path.append("../bin/")
import pyHiChi as hichi
import numpy as np
from numba import cfunc

# -------------- fields and currents --------------
# all functions are precompiled for the best performance

@cfunc("float64(float64,float64,float64)")  
def null_value(x, y, z):  # start conditions
    return 0.0
    
@cfunc("float64(float64,float64,float64)")
def block(x, a, b):
    return 1.0 if a <= x < b else 0.0
    
@cfunc("float64(float64,float64,float64,float64)")  
def func_jx(x, y, z, t):
    return 0.0
    
@cfunc("float64(float64,float64,float64,float64)")  
def func_jy(x, y, z, t):
    return 0.0
    
@cfunc("float64(float64,float64,float64,float64)")  
def func_jz(x, y, z, t):
    T = 8
    Tc = 4 * hichi.c
    return np.sin(2.0*np.pi*t/T) * \
        np.cos(np.pi*x/Tc)**2 * block(x, -0.5*Tc, 0.5*Tc) * \
        np.cos(np.pi*y/Tc)**2 * block(y, -0.5*Tc, 0.5*Tc)
    
# -------------- simulation parameters --------------

grid_size = hichi.Vector3d(128, 128, 1)
grid_step = hichi.Vector3d(hichi.c, hichi.c, hichi.c)

min_coords = hichi.Vector3d(-grid_size.x/2*hichi.c, -grid_size.y/2*hichi.c, 0.0)
max_coords = hichi.Vector3d(grid_size.x/2*hichi.c, grid_size.y/2*hichi.c, 0.0)

time_step = 0.4*grid_step.x / hichi.c

field = hichi.PSATDField(grid_size, min_coords, grid_step, time_step)
#field = hichi.YeeField(grid_size, min_coords, grid_step, time_step)
field.set_E(null_value.address, null_value.address, null_value.address)
field.set_B(null_value.address, null_value.address, null_value.address)

def update_fields(t):
    field.set_J(func_jx.address, func_jy.address, func_jz.address, t)
    field.update_fields()

# -------------- show --------------

import matplotlib.pyplot as plt
import matplotlib.animation as animation

N = 128
x = np.linspace(min_coords.x, max_coords.x, N)
y = np.linspace(min_coords.y, max_coords.y, N)

def get_fields():
    global field, x, y, N
    Ez = np.zeros(shape=(N,N))
    for ix in range(N):
        for iy in range(N):
            coord_xy = hichi.Vector3d(x[ix], y[iy], 0.0)
            Ez[ix, iy] = field.get_Ez(coord_xy)
    return Ez 

Ez = get_fields()

fig, ax = plt.subplots(ncols=1, nrows=1)

im = ax.imshow(Ez, cmap='RdBu', interpolation='none',
               extent=(min_coords.x, max_coords.x, min_coords.y, max_coords.y),
               animated=True)
fig.colorbar(im, ax=ax)
ax.set_title("Ez")
ax.set_xlabel("x")
ax.set_ylabel("y")

def update_fig(iter):
    t = iter * time_step
    update_fields(t)
    Ez = get_fields()
    im.set_array(Ez)
    return im, 
    
ani = animation.FuncAnimation(fig, update_fig, interval=10, blit=True)

plt.tight_layout()

plt.show()
