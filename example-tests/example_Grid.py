import sys
sys.path.append("../bin/")
import pyHiChi as pfc
import numpy as np
import math as ma

def valueE(x, y, z):
    E = pfc.vector3d(0, np.cos(z), 0) #sin(x)
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
    B = pfc.vector3d(-np.cos(z), 0, 0)
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

def step(minCoords, maxCoords, gridSize):
    steps = pfc.vector3d(1, 1, 1)
    steps.x = (maxCoords.x - minCoords.x)/(gridSize.x)
    steps.y = (maxCoords.y - minCoords.y)/(gridSize.y)
    steps.z = (maxCoords.z - minCoords.z)/(gridSize.z)
    return steps

gridSize = pfc.vector3d(5, 10, 11)
minCoords = pfc.vector3d(0.0, 1.0, 0.0)
maxCoords = pfc.vector3d(3.5, 7.0, 2*ma.pi)
stepsGrid = step(minCoords, maxCoords, gridSize)
timeStep = 1e-16

grid1 = pfc.YeeGrid(gridSize, timeStep, minCoords, stepsGrid)
grid2 = pfc.YeeGrid(gridSize, timeStep, minCoords, stepsGrid)
    
grid1.setE(valueE)
grid1.setB(valueB)

grid2.setE(valueEx, valueEy, valueEz)
grid2.setB(valueBx, valueBy, valueBz)


#show
import matplotlib.pyplot as plt

N = 37
x = np.arange(0, 3.5, 3.5/N)
z = np.arange(0, 2*ma.pi, 2*ma.pi/N)

Ex1 = np.zeros(shape=(N,N))
Ex2 = np.zeros(shape=(N,N))
Ey1 = np.zeros(shape=(N,N))
Ey2 = np.zeros(shape=(N,N))
Bx1 = np.zeros(shape=(N,N))
Bx2 = np.zeros(shape=(N,N))

for ix in range(N):
    for iy in range(N):
        coordXZ = pfc.vector3d(x[ix], 0.0, z[iy])
        E1 = grid1.getE(coordXZ)
        Ex1[ix, iy] = E1.x
        Ey1[ix, iy] = E1.y
        Bx1[ix, iy] = grid1.getB(coordXZ).x
        E2 = grid2.getE(coordXZ)
        Ex2[ix, iy] = E2.x
        Ey2[ix, iy] = E2.y
        Bx2[ix, iy] = grid2.getB(coordXZ).x


fig, axes = plt.subplots(ncols=3, nrows=2)

bar11 = axes[0, 0].imshow(Ex1, cmap='RdBu', interpolation='none', extent=(0, 2*ma.pi, 0, 3.5))
fig.colorbar(bar11, ax=axes[0, 0])
axes[0, 0].set_title("Ex1")
axes[0, 0].set_xlabel("x")
axes[0, 0].set_ylabel("z")

bar12 = axes[0, 1].imshow(Ey1, cmap='RdBu', interpolation='none', extent=(0, 2*ma.pi, 0, 3.5))
fig.colorbar(bar12, ax=axes[0, 1])
axes[0, 1].set_title("Ey1")
axes[0, 1].set_xlabel("x")
axes[0, 1].set_ylabel("z")

bar13 = axes[0, 2].imshow(Bx1, cmap='RdBu', interpolation='none', extent=(0, 2*ma.pi, 0, 3.5))
fig.colorbar(bar13, ax=axes[0, 2])
axes[0, 2].set_title("Bx1")
axes[0, 2].set_xlabel("x")
axes[0, 2].set_ylabel("z")

bar21 = axes[1, 0].imshow(Ex2, cmap='RdBu', interpolation='none', extent=(0, 2*ma.pi, 0, 3.5))
fig.colorbar(bar21, ax=axes[1, 0])
axes[1, 0].set_title("Ex2")
axes[1, 0].set_xlabel("x")
axes[1, 0].set_ylabel("z")

bar22 = axes[1, 1].imshow(Ey2, cmap='RdBu', interpolation='none', extent=(0, 2*ma.pi, 0, 3.5))
fig.colorbar(bar22, ax=axes[1, 1])
axes[1, 1].set_title("Ey2")
axes[1, 1].set_xlabel("x")
axes[1, 1].set_ylabel("z")

bar23 = axes[1, 2].imshow(Bx2, cmap='RdBu', interpolation='none', extent=(0, 2*ma.pi, 0, 3.5))
cbar = fig.colorbar(bar23, ax=axes[1, 2])
axes[1, 2].set_title("Bx2")
axes[1, 2].set_xlabel("x")
axes[1, 2].set_ylabel("z")

plt.tight_layout()

plt.show()


