import pyHiChi as hichi
import numpy as np
import math as ma
import matplotlib.pyplot as plt

def getFields(grid, minCoords, maxCoords):
    N = 128
    x = np.linspace(minCoords.x, maxCoords.x, N)
    y = np.linspace(minCoords.y, maxCoords.y, N)
    Ey = np.zeros(shape=(N,N))
    for ix in range(N):
        for iy in range(N):
            coordXY = hichi.vector3d(x[ix], y[iy], 0.0)
            E = grid.getE(coordXY)
            Ey[N - iy - 1, ix] = E.y
    return Ey


# ---------- example1 ------------------
# set -> grid -> inverse mapping -> get
# next code makes rotation of a pulse

def fieldValue1(x, y, z):
    return np.exp(-x**2-y**2-z**2)*np.sin(5*x)
    
def nullValue(x, y, z):
    return 0.0
    
minCoords = hichi.vector3d(-5, -5, -5)
maxCoords = hichi.vector3d(5, 5, 5)

rotationZMapping = hichi.RotationMapping(hichi.Axis.z, ma.pi/3)

gridSize = hichi.vector3d(64, 64, 2)

gridStep = hichi.vector3d(
            (maxCoords.x - minCoords.x) / gridSize.x,
            (maxCoords.y - minCoords.y) / gridSize.y,
            (maxCoords.z - minCoords.z) / gridSize.z
           )
           
timeStep = 0.0

# create the computational grid which keeps not transformed fields
grid1 = hichi.PSATDGrid(gridSize, timeStep, minCoords, gridStep)

# create shallow copy of the grid which makes the transform when getter is called
gridMapping1 = hichi.PSATDGridMapping(grid1)
gridMapping1.setMapping(rotationZMapping)

# set direct fields to the computational grid
grid1.setE(nullValue, fieldValue1, nullValue)
grid1.setB(nullValue, nullValue, fieldValue1)

Ey1Direct = getFields(grid1, minCoords, maxCoords)  # this is not transformed field
Ey1Transformed = getFields(gridMapping1, minCoords, maxCoords)  # this is trasformed field


# ---------- example2 -----------------------------------
# set -> direct mapping -> grid -> inverse mapping -> get
# next code uses transform that is used for tight focusing

R0 = 3.5
L = 1.0

# just a piece of sphere
def fieldValue2(x, y, z):

    def lf(x):
        return 1.0 if x >= -L*0.5 and x <= L*0.5 else 0.0
    
    def tsf(angle):
        return 1.0 if angle <= 0.35*ma.pi else 0.0   
      
    R = np.sqrt(x*x + y*y + z*z)
    if (R > 1e-5):
        angle = np.arcsin(np.sqrt(y*y + z*z)/R)
        return (1.0/R) * lf(R - R0) * tsf(angle) * (x < 0.0)
    else:
        return 0.0     

D = 2  # width of the band along x axis that will be computed

mappingTF = hichi.TightFocusingMapping(R0, L, D)

xMin = mappingTF.getxMin()
xMax = mappingTF.getxMax()  # xMax = xMin + D

minInverseCoords = hichi.vector3d(xMin, minCoords.y, minCoords.z)
maxInverseCoords = hichi.vector3d(xMax, maxCoords.y, maxCoords.z)
gridStep = hichi.vector3d((maxInverseCoords.x - minInverseCoords.x) / gridSize.x, \
                          (maxInverseCoords.y - minInverseCoords.y) / gridSize.y, \
                          (maxInverseCoords.z - minInverseCoords.z) / gridSize.z)

# create the computational grid
grid2 = hichi.PSATDGrid(gridSize, timeStep, minInverseCoords, gridStep)

# create shallow copy of the grid
gridMapping2 = hichi.PSATDGridMapping(grid2)
gridMapping2.setMapping(mappingTF)

# lines 96 and 99 can be replaced by only one line
# gridMapping2 = hichi.PSATDGridMapping(gridSize, timeStep, minInverseCoords, gridStep)
# so it will not possible to get fields that are keeped in the computational grid

# set fields to the grid with mapping, so there will make direct transform before setting
gridMapping2.setE(nullValue, fieldValue2, nullValue)
gridMapping2.setB(nullValue, nullValue, fieldValue2)

Ey2Transformed = getFields(grid2, minInverseCoords, maxInverseCoords)  # this field is keeped in the computational grid
Ey2Direct = getFields(gridMapping2, minCoords, maxCoords)  # this field is getted by making inverse transform for fields in grid
                                                           # so we will get "corners" near the y axis


# -------------------- show ----------------------------


fig, axes = plt.subplots(ncols=2, nrows=2)


im1 = axes[0, 0].imshow(Ey1Direct, cmap='RdBu', interpolation='none',
    extent=(minCoords.x, maxCoords.x, minCoords.y, maxCoords.y), animated = True)
fig.colorbar(im1, ax=axes[0, 0])
axes[0, 0].set_xlabel("x")
axes[0, 0].set_ylabel("y")
axes[0, 0].set_title("Ey1Direct")


im2 = axes[0, 1].imshow(Ey1Transformed, cmap='RdBu', interpolation='none',
    extent=(minCoords.x, maxCoords.x, minCoords.y, maxCoords.y), animated = True)
fig.colorbar(im2, ax=axes[0, 1])
axes[0, 1].set_xlabel("x")
axes[0, 1].set_ylabel("y")
axes[0, 1].set_title("Ey1Transformed")


im3 = axes[1, 0].imshow(Ey2Direct, cmap='RdBu', interpolation='none',
    extent=(minCoords.x, maxCoords.x, minCoords.y, maxCoords.y), animated = True)
fig.colorbar(im3, ax=axes[1, 0])
axes[1, 0].set_xlabel("x")
axes[1, 0].set_ylabel("y")
axes[1, 0].set_title("Ey2Direct")


im4 = axes[1, 1].imshow(Ey2Transformed, cmap='RdBu', interpolation='none',
    extent=(minInverseCoords.x, maxInverseCoords.x, minInverseCoords.y, maxInverseCoords.y), animated = True)
fig.colorbar(im4, ax=axes[1, 1])
axes[1, 1].set_xlabel("x")
axes[1, 1].set_ylabel("y")
axes[1, 1].set_title("Ey2Transformed")


fig.tight_layout()

plt.show()
