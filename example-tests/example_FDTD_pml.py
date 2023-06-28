import sys
sys.path.append("../bin/")
import pyHiChi as hichi
import numpy as np
import math as ma

min_coord = 0.0
max_coord = 1.0

def block(x, x_min, x_max):
    return 1.0 if x_min <= x < x_max else 0.0

def valueEx(x, y, z):
    Ex = np.sin((z - min_coord) * 2.0 * hichi.pi / (max_coord - min_coord)) * \
        block(z, min_coord, max_coord)
    return Ex

def valueEy(x, y, z):
    Ey = 0
    return Ey

def valueEz(x, y, z):
    Ez = 0
    return Ez

def valueBx(x, y, z):
    Bx = 0
    return Bx

def valueBy(x, y, z):
    By = -np.sin((z - min_coord) * 2.0 * hichi.pi / (max_coord - min_coord)) * \
        block(z, min_coord, max_coord)
    return By

def valueBz(x, y, z):
    Bz = 0
    return Bz

grid_size = hichi.Vector3d(16, 16, 64)
d_coord = (max_coord - min_coord) / grid_size.z
grid_step = hichi.Vector3d(d_coord, d_coord, d_coord)
pml_size = hichi.Vector3d(0, 0, 8)

min_coords = hichi.Vector3d(0.0, 0.0, 0.0) - pml_size*grid_step
max_coords = min_coords + grid_step*(grid_size + 2*pml_size)

time_step = 0.1*grid_step.norm() / hichi.c

field = hichi.YeeField(grid_size + 2*pml_size, min_coords, grid_step, time_step)
field.set_E(valueEx, valueEy, valueEz)
field.set_B(valueBx, valueBy, valueBz)

field.set_PML(pml_size)
field.set_periodical_BC(hichi.Axis.X)
field.set_periodical_BC(hichi.Axis.Y)

#show
import matplotlib.pyplot as plt
import matplotlib.animation as animation

N = 100
x = np.linspace(min_coords.x, max_coords.x, N)
z = np.linspace(min_coords.z, max_coords.z, N)


def get_fields():
    global field, x, z, N
    #print(field)
    Ex = np.zeros(shape=(N,N))
    Ey = np.zeros(shape=(N,N))
    Ez = np.zeros(shape=(N,N))
    Bx = np.zeros(shape=(N,N))
    By = np.zeros(shape=(N,N))
    Bz = np.zeros(shape=(N,N))
    for ix in range(N):
        for iy in range(N):
            coord_xz = hichi.Vector3d(x[ix], 0.0, z[iy])
            E = field.get_E(coord_xz)
            Ex[ix, iy] = E.x
            Ey[ix, iy] = E.y
            Ez[ix, iy] = E.z
            B = field.get_B(coord_xz)
            Bx[ix, iy] = B.x
            By[ix, iy] = B.y
            Bz[ix, iy] = B.z
    return Ex, Ey, Ez, Bx, By, Bz

def update_data():
    for i in range(10):
        field.update_fields()  

(Ex, Ey, Ez, Bx, By, Bz) = get_fields()
fig, axes = plt.subplots(ncols=3, nrows=2)

im11 = axes[0, 0].imshow(Ex, cmap='RdBu', interpolation='none',
    extent=(min_coords.z, max_coords.z, min_coords.x, max_coords.x), animated = True)
fig.colorbar(im11, ax=axes[0, 0])
axes[0, 0].set_title("Ex")
axes[0, 0].set_xlabel("z")
axes[0, 0].set_ylabel("x")

im12 = axes[0, 1].imshow(Ey, cmap='RdBu', interpolation='none',
    extent=(min_coords.z, max_coords.z, min_coords.x, max_coords.x), animated = True)
fig.colorbar(im12, ax=axes[0, 1])
axes[0, 1].set_title("Ey")
axes[0, 1].set_xlabel("z")
axes[0, 1].set_ylabel("x")

im13 = axes[0, 2].imshow(Ez, cmap='RdBu', interpolation='none',
    extent=(min_coords.z, max_coords.z, min_coords.x, max_coords.x), animated = True)
fig.colorbar(im13, ax=axes[0, 2])
axes[0, 2].set_title("Ez")
axes[0, 2].set_xlabel("z")
axes[0, 2].set_ylabel("x")

im21 = axes[1, 0].imshow(Bx, cmap='RdBu', interpolation='none',
    extent=(min_coords.z, max_coords.z, min_coords.x, max_coords.x), animated = True)
fig.colorbar(im21, ax=axes[1, 0])
axes[1, 0].set_title("Bx")
axes[1, 0].set_xlabel("z")
axes[1, 0].set_ylabel("x")

im22 = axes[1, 1].imshow(By, cmap='RdBu', interpolation='none',
    extent=(min_coords.z, max_coords.z, min_coords.x, max_coords.x), animated = True)
fig.colorbar(im22, ax=axes[1, 1])
axes[1, 1].set_title("By")
axes[1, 1].set_xlabel("z")
axes[1, 1].set_ylabel("x")

im23 = axes[1, 2].imshow(Bz, cmap='RdBu', interpolation='none',
    extent=(min_coords.z, max_coords.z, min_coords.x, max_coords.x), animated = True)
fig.colorbar(im23, ax=axes[1, 2])
axes[1, 2].set_title("Bz")
axes[1, 2].set_xlabel("z")
axes[1, 2].set_ylabel("x")

def compute_energy(Ex, Ey, Ez, Bx, By, Bz):
    return np.sum(Ex**2 + Ey**2 + Ez**2 + Bx**2 + By**2 + Bz**2)

def update_fig(*args):
    update_data()
    (Ex, Ey, Ez, Bx, By, Bz) = get_fields()
    im11.set_array(Ex)
    im12.set_array(Ey)
    im13.set_array(Ez)
    im21.set_array(Bx)
    im22.set_array(By)
    im23.set_array(Bz)
    print(compute_energy(Ex, Ey, Ez, Bx, By, Bz))
    return im11, im12, im13, im21, im22, im23, 
    
ani = animation.FuncAnimation(fig, update_fig, interval=50, blit=True)

plt.tight_layout()

plt.show()


