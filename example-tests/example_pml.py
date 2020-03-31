import pyHiChi as pfc
import numpy as np
import math as ma

gridSize = pfc.vector3d(144, 30, 1)
pmlSize = pfc.vector3d(8, 0, 0)
minCoords = pfc.vector3d(0.0, 0.0, 0.0)
maxCoords = pfc.vector3d(gridSize.x * pfc.c, gridSize.y * pfc.c, gridSize.z * pfc.c)

def step(minCoords, maxCoords, gridSize):
    steps = pfc.vector3d(1, 1, 1)
    steps.x = (maxCoords.x - minCoords.x)/(gridSize.x)
    steps.y = (maxCoords.y - minCoords.y)/(gridSize.y)
    steps.z = (maxCoords.z - minCoords.z)/(gridSize.z)
    return steps
    
stepsGrid = step(minCoords, maxCoords, gridSize)
timeStep = 0.1

pmlLeftEnd = minCoords.x+pmlSize.x*stepsGrid.x
pmlRightStart = maxCoords.x-pmlSize.x*stepsGrid.x
internalWidth = pmlRightStart-pmlLeftEnd#maxCoords.x-minCoords.x

def valueEx(x, y, z):
    Ex = 0
    return Ex
def valueEy(x, y, z):
    if (x<pmlLeftEnd or x>=pmlRightStart):
        Ey=0
    else: 
        Ey = ma.sin(2*ma.pi/internalWidth*(x-pmlLeftEnd))
    return Ey
def valueEz(x, y, z):
    Ez = 0
    return Ez

def valueBx(x, y, z):
    Bx = 0
    return Bx
def valueBy(x, y, z):
    By = 0
    return By
def valueBz(x, y, z):
    if (x<pmlLeftEnd or x>=pmlRightStart):
        Bz=0
    else: 
        Bz = ma.sin(2*ma.pi/internalWidth*(x-pmlLeftEnd))
    return Bz

grid = pfc.PSTDGrid(gridSize, timeStep, minCoords, stepsGrid)
# grid = pfc.PSATDTimeStraggeredGrid(gridSize, timeStep, minCoords, stepsGrid)
# grid = pfc.YeeGrid(gridSize, timeStep, minCoords, stepsGrid) 
grid.setE(valueEx, valueEy, valueEz)
grid.setB(valueBx, valueBy, valueBz)

fieldSolver = pfc.PSTD(grid)
# fieldSolver = pfc.PSATDTimeStraggered(grid)
# fieldSolver = pfc.FDTD(grid)
fieldSolver.setPML(int(pmlSize.x), int(pmlSize.y), int(pmlSize.z))

#show
import matplotlib.pyplot as plt
import matplotlib.animation as animation

Nx = int(gridSize.x)
Ny = int(gridSize.y)
x = np.arange(minCoords.x, maxCoords.x, (maxCoords.x-minCoords.x)/Nx)
y = np.arange(minCoords.y, maxCoords.y, (maxCoords.y-minCoords.y)/Ny)

def getFields():
    global grid, x, y, Nx, Ny
    Ey = np.zeros(shape=(Ny,Nx))
    Bz = np.zeros(shape=(Ny,Nx))
    for ix in range(Nx):
        for iy in range(Ny):
            coordXY = pfc.vector3d(x[ix], y[iy], 0.0)
            E = grid.getE(coordXY)
            Ey[iy, ix] = E.y
            B = grid.getB(coordXY)
            Bz[iy, ix] = B.z
    return Ey, Bz

def updateData():
    for i in range(3):
        fieldSolver.updateFields()
        
def computeEnergy():
    (Ey, Bz) = getFields()
    energy=0
    for ix in range(Nx):
        for iy in range(Ny):
            energy+=Ey[iy, ix]**2+Bz[iy, ix]**2
    return energy

(Ey, Bz) = getFields()

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

iter=0
def updatefig(*args):
    global iter
    updateData()
    (Ey, Bz) = getFields()
    im11.set_array(Ey)
    im12.set_array(Bz)
    if (iter%50 == 0):
        print("Energy = "+str(computeEnergy()))
    iter+=1
    return im11, im12
    
ani = animation.FuncAnimation(fig, updatefig, interval=10, blit=True)

plt.tight_layout()

plt.show()





