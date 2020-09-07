import sys
sys.path.append("../bin/")
import pyHiChi as hichi
import numpy as np

def valueE(x, y, z):
    E = hichi.Vector3d(0, np.cos(z), 0) #sin(x)
    return E

def valueEx(x, y, z):
    Ex = 0
    return Ex

def valueEy(x, y, z):
    Ey = np.cos(z)
    return Ey

def valueEz(x, y, z):
    Ez = 0
    return Ez


def valueB(x, y, z):
    B = hichi.Vector3d(-np.cos(z), 0, 0)
    return B

def valueBx(x, y, z):
    Bx = -np.cos(z)
    return Bx

def valueBy(x, y, z):
    By = 0
    return By

def valueBz(x, y, z):
    Bz = 0
    return Bz


field_size = hichi.Vector3d(5, 10, 11)
min_coords = hichi.Vector3d(0.0, 1.0, 0.0)
max_coords = hichi.Vector3d(3.5, 7.0, 2*np.pi)
field_step = (max_coords - min_coords) / field_size
time_step = 1e-16

field1 = hichi.YeeField(field_size, min_coords, field_step, time_step)
field2 = hichi.YeeField(field_size, min_coords, field_step, time_step)
    
field1.set_E(valueE)
field1.set_B(valueB)

field2.set_E(valueEx, valueEy, valueEz)
field2.set_B(valueBx, valueBy, valueBz)


#show
import matplotlib.pyplot as plt

N = 37
x = np.arange(0, 3.5, 3.5/N)
z = np.arange(0, 2*np.pi, 2*np.pi/N)

Ex1 = np.zeros(shape=(N,N))
Ex2 = np.zeros(shape=(N,N))
Ey1 = np.zeros(shape=(N,N))
Ey2 = np.zeros(shape=(N,N))
Bx1 = np.zeros(shape=(N,N))
Bx2 = np.zeros(shape=(N,N))

for ix in range(N):
    for iy in range(N):
        coord_xz = hichi.Vector3d(x[ix], 0.0, z[iy])
        E1 = field1.get_E(coord_xz)
        Ex1[ix, iy] = E1.x
        Ey1[ix, iy] = E1.y
        Bx1[ix, iy] = field1.get_B(coord_xz).x
        E2 = field2.get_E(coord_xz)
        Ex2[ix, iy] = E2.x
        Ey2[ix, iy] = E2.y
        Bx2[ix, iy] = field2.get_B(coord_xz).x


fig, axes = plt.subplots(ncols=3, nrows=2)

bar11 = axes[0, 0].imshow(Ex1, cmap='RdBu', interpolation='none', extent=(0, 2*np.pi, 0, 3.5))
fig.colorbar(bar11, ax=axes[0, 0])
axes[0, 0].set_title("Ex1")
axes[0, 0].set_xlabel("x")
axes[0, 0].set_ylabel("z")

bar12 = axes[0, 1].imshow(Ey1, cmap='RdBu', interpolation='none', extent=(0, 2*np.pi, 0, 3.5))
fig.colorbar(bar12, ax=axes[0, 1])
axes[0, 1].set_title("Ey1")
axes[0, 1].set_xlabel("x")
axes[0, 1].set_ylabel("z")

bar13 = axes[0, 2].imshow(Bx1, cmap='RdBu', interpolation='none', extent=(0, 2*np.pi, 0, 3.5))
fig.colorbar(bar13, ax=axes[0, 2])
axes[0, 2].set_title("Bx1")
axes[0, 2].set_xlabel("x")
axes[0, 2].set_ylabel("z")

bar21 = axes[1, 0].imshow(Ex2, cmap='RdBu', interpolation='none', extent=(0, 2*np.pi, 0, 3.5))
fig.colorbar(bar21, ax=axes[1, 0])
axes[1, 0].set_title("Ex2")
axes[1, 0].set_xlabel("x")
axes[1, 0].set_ylabel("z")

bar22 = axes[1, 1].imshow(Ey2, cmap='RdBu', interpolation='none', extent=(0, 2*np.pi, 0, 3.5))
fig.colorbar(bar22, ax=axes[1, 1])
axes[1, 1].set_title("Ey2")
axes[1, 1].set_xlabel("x")
axes[1, 1].set_ylabel("z")

bar23 = axes[1, 2].imshow(Bx2, cmap='RdBu', interpolation='none', extent=(0, 2*np.pi, 0, 3.5))
cbar = fig.colorbar(bar23, ax=axes[1, 2])
axes[1, 2].set_title("Bx2")
axes[1, 2].set_xlabel("x")
axes[1, 2].set_ylabel("z")

plt.tight_layout()

plt.show()


