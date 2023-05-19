import sys
sys.path.append("../bin/")
import pyHiChi as hichi
import numpy as np
from numba import cfunc

# -------------- fields and currents --------------
# all functions are precompiled for better performance

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

grid_size = hichi.Vector3d(256, 256, 1)
grid_step = hichi.Vector3d(hichi.c, hichi.c, hichi.c)

min_coords = hichi.Vector3d(-grid_size.x/2*hichi.c, -grid_size.y/2*hichi.c, 0.0)
max_coords = hichi.Vector3d(grid_size.x/2*hichi.c, grid_size.y/2*hichi.c, 0.0)

time_step = 0.4*grid_step.x / hichi.c

field = hichi.PSATDField(grid_size, min_coords, grid_step, time_step)
# field = hichi.YeeField(grid_size, min_coords, grid_step, time_step)
field.set_E(null_value.address, null_value.address, null_value.address)
field.set_B(null_value.address, null_value.address, null_value.address)
# field.set_PML(8, 8, 0)  # set pml if necessary

def update_fields(t):
    field.set_J(func_jx.address, func_jy.address, func_jz.address, t)
    field.update_fields()


# -------------- show --------------

import matplotlib.pyplot as plt
import matplotlib.animation as animation

N = 512
x = np.linspace(min_coords.x, max_coords.x, N)
y = np.linspace(min_coords.y, max_coords.y, N)
    
def get_fields_1d():
    return field.get_Ez_x_line(z_pos=0.0, y_pos=0.0,
        x_min=min_coords.x, x_max=max_coords.x, x_size=N,)

def get_fields_2d():
    return field.get_Ez_xy_plane(z_pos=0.0,
        x_min=min_coords.x, x_max=max_coords.x, x_size=N,
        y_min=min_coords.y, y_max=max_coords.y, y_size=N,)

fig, ax = plt.subplots(ncols=2, nrows=1)

im = ax[0].imshow(get_fields_2d(), cmap='RdBu', interpolation='none',
                  extent=(min_coords.x, max_coords.x, min_coords.y, max_coords.y),
                  animated=True, vmin=-3, vmax=3)
fig.colorbar(im, ax=ax[0])
ax[0].set_title("Ez")
ax[0].set_xlabel("x")
ax[0].set_ylabel("y")

line = ax[1].plot(x, get_fields_1d())[0]
ax[1].set_title("Ez")
ax[1].set_xlabel("x")
ax[1].grid()
ax[1].set_ylim((-10, 10))

def update_fig(iter):
    t = iter * time_step
    update_fields(t)
    im.set_array(get_fields_2d())
    line.set_data(x, get_fields_1d())
    return im, line,
    
ani = animation.FuncAnimation(fig, update_fig, interval=10, blit=True)

plt.tight_layout()

plt.show()
