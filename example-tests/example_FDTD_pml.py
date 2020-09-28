import sys
sys.path.append("../bin/")
import pyHiChi as hichi
import numpy as np
import math as ma

def valueEx(x, y, z):
    Ex = 0
    return Ex

def valueEy(x, y, z):
    Ey = np.cos(z + np.pi/6)
    return Ey

def valueEz(x, y, z):
    Ez = 0
    return Ez


def valueBx(x, y, z):
    Bx = -np.cos(z + np.pi/6)
    return Bx

def valueBy(x, y, z):
    By = 0
    return By

def valueBz(x, y, z):
    Bz = 0
    return Bz

grid_size = hichi.Vector3d(40, 40, 40)
min_coords = hichi.Vector3d(0.0, 1.0, 0.0)
max_coords = hichi.Vector3d(3.5, 7.0, 2*np.pi)
grid_step = (max_coords - min_coords) / grid_size
time_step = 1e-13

field = hichi.YeeField(grid_size, min_coords, grid_step, time_step)
field.set_E(valueEx, valueEy, valueEz)
field.set_B(valueBx, valueBy, valueBz)

field.set_PML(0, 0, 0)

#show
import matplotlib.pyplot as plt
import matplotlib.animation as animation

N = 30
x = np.arange(0, 3.5, 3.5/N)
z = np.arange(0, 2*np.pi, 2*np.pi/N)


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

print(Ey[10, 10])
fig, axes = plt.subplots(ncols=3, nrows=2)

im11 = axes[0, 0].imshow(Ex, cmap='RdBu', interpolation='none', extent=(0, 2*np.pi, 0, 3.5), animated = True)
fig.colorbar(im11, ax=axes[0, 0])
axes[0, 0].set_title("Ex")
axes[0, 0].set_xlabel("x")
axes[0, 0].set_ylabel("z")

im12 = axes[0, 1].imshow(Ey, cmap='RdBu', interpolation='none', extent=(0, 2*np.pi, 0, 3.5), animated = True)
fig.colorbar(im12, ax=axes[0, 1])
axes[0, 1].set_title("Ey")
axes[0, 1].set_xlabel("x")
axes[0, 1].set_ylabel("z")

im13 = axes[0, 2].imshow(Ez, cmap='RdBu', interpolation='none', extent=(0, 2*np.pi, 0, 3.5), animated = True)
fig.colorbar(im13, ax=axes[0, 2])
axes[0, 2].set_title("Ez")
axes[0, 2].set_xlabel("x")
axes[0, 2].set_ylabel("z")

im21 = axes[1, 0].imshow(Bx, cmap='RdBu', interpolation='none', extent=(0, 2*np.pi, 0, 3.5), animated = True)
fig.colorbar(im21, ax=axes[1, 0])
axes[1, 0].set_title("Bx")
axes[1, 0].set_xlabel("x")
axes[1, 0].set_ylabel("z")

im22 = axes[1, 1].imshow(By, cmap='RdBu', interpolation='none', extent=(0, 2*np.pi, 0, 3.5), animated = True)
fig.colorbar(im22, ax=axes[1, 1])
axes[1, 1].set_title("By")
axes[1, 1].set_xlabel("x")
axes[1, 1].set_ylabel("z")

im23 = axes[1, 2].imshow(Bz, cmap='RdBu', interpolation='none', extent=(0, 2*np.pi, 0, 3.5), animated = True)
fig.colorbar(im23, ax=axes[1, 2])
axes[1, 2].set_title("Bz")
axes[1, 2].set_xlabel("x")
axes[1, 2].set_ylabel("z")


def update_fig(*args):
    update_data()
    (Ex, Ey, Ez, Bx, By, Bz) = get_fields()
    im11.set_array(Ex)
    im12.set_array(Ey)
    im13.set_array(Ez)
    im21.set_array(Bx)
    im22.set_array(By)
    im23.set_array(Bz)
    print(Ey[10, 10])
    return im11, im12, im13, im21, im22, im23, 
    
ani = animation.FuncAnimation(fig, update_fig, interval=50, blit=True)

plt.tight_layout()

plt.show()


