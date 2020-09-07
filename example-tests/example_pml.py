import sys
sys.path.append("../bin/")
import pyHiChi as hichi
import numpy as np

field_size = hichi.Vector3d(144, 30, 1)
pml_size = hichi.Vector3d(8, 0, 0)
min_coords = hichi.Vector3d(0.0, 0.0, 0.0)
max_coords = hichi.Vector3d(field_size.x * hichi.c, field_size.y * hichi.c, field_size.z * hichi.c)
    
field_step = (max_coords - min_coords) / field_size
time_step = 0.1

pml_left_end = min_coords.x + pml_size.x*field_step.x
pml_right_end = max_coords.x - pml_size.x*field_step.x
internal_width = pml_right_end - pml_left_end

def value_Ex(x, y, z):
    Ex = 0
    return Ex

def value_Ey(x, y, z):
    if (x < pml_left_end or x >= pml_right_end):
        Ey=0
    else: 
        Ey = np.sin(2*np.pi/internal_width*(x - pml_left_end))
    return Ey

def value_Ez(x, y, z):
    Ez = 0
    return Ez


def value_Bx(x, y, z):
    Bx = 0
    return Bx

def value_By(x, y, z):
    By = 0
    return By

def value_Bz(x, y, z):
    if (x < pml_left_end or x >= pml_right_end):
        Bz=0
    else: 
        Bz = np.sin(2*np.pi/internal_width*(x - pml_left_end))
    return Bz

field = hichi.PSTDField(field_size, min_coords, field_step, time_step)
# field = hichi.PSATDField(field_size, min_coords, field_step, time_step)
# field = hichi.PSATDTimeStraggeredField(field_size, min_coords, field_step, time_step)
# field = hichi.YeeField(field_size, min_coords, field_step, time_step) 
field.set_E(value_Ex, value_Ey, value_Ez)
field.set_B(value_Bx, value_By, value_Bz)
field.set_PML(int(pml_size.x), int(pml_size.y), int(pml_size.z))


# ------------------- show -------------------

import matplotlib.pyplot as plt
import matplotlib.animation as animation

Nx = int(field_size.x)
Ny = int(field_size.y)
x = np.arange(min_coords.x, max_coords.x, (max_coords.x - min_coords.x)/Nx)
y = np.arange(min_coords.y, max_coords.y, (max_coords.y - min_coords.y)/Ny)

def get_fields():
    global field, x, y, Nx, Ny
    Ey = np.zeros(shape=(Ny, Nx))
    Bz = np.zeros(shape=(Ny, Nx))
    for ix in range(Nx):
        for iy in range(Ny):
            coord_xy = hichi.Vector3d(x[ix], y[iy], 0.0)
            E = field.get_E(coord_xy)
            Ey[iy, ix] = E.y
            B = field.get_B(coord_xy)
            Bz[iy, ix] = B.z
    return Ey, Bz

def update_data():
    for i in range(3):
        field.update_fields()
        
def compute_energy():
    (Ey, Bz) = get_fields()
    energy=0
    for ix in range(Nx):
        for iy in range(Ny):
            energy += Ey[iy, ix]**2 + Bz[iy, ix]**2
    return energy

(Ey, Bz) = get_fields()

fig, axes = plt.subplots(ncols=2, nrows=1)

im11 = axes[0].imshow(Ey, cmap='RdBu', interpolation='none', extent=(0, 2, 0, 1), animated = True)
fig.colorbar(im11, ax=axes[0])
axes[0].set_title("Ey")
axes[0].set_xlabel("x")
axes[0].set_ylabel("y")

im12 = axes[1].imshow(Bz, cmap='RdBu', interpolation='none', extent=(0, 2, 0, 1), animated = True)
fig.colorbar(im12, ax=axes[1])
axes[1].set_title("Bz")
axes[1].set_xlabel("x")
axes[1].set_ylabel("y")

iter = 0
def update_fig(*args):
    global iter
    update_data()
    (Ey, Bz) = get_fields()
    im11.set_array(Ey)
    im12.set_array(Bz)
    if (iter % 50 == 0):
        print("Energy = " + str(compute_energy()))
    iter += 1
    return im11, im12
    
ani = animation.FuncAnimation(fig, update_fig, interval=10, blit=True)

plt.tight_layout()

plt.show()





