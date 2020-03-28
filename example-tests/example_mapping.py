import pyHiChi as hichi
import numpy as np
import math as ma


def fieldValue(x, y, z):
    return np.exp(-x**2-y**2-z**2)*np.sin(5*x)
    
def nullValue(x, y, z):
    return 0.0
    
minDirectCoords = hichi.vector3d(-5, -5, -5)
maxDirectCoords = hichi.vector3d(5, 5, 5)

# create some transforms
angle = ma.pi/3
rotationZMapping = hichi.RotationMapping(hichi.Axis.z, angle)

angle2 = ma.pi/6
rotationZMapping2 = hichi.RotationMapping(hichi.Axis.z, angle2)

shift = hichi.vector3d(-2.0, -2.0, 0.0)
shiftMapping = hichi.ShiftMapping(shift)

scaleCoef = 1.5
scaleMapping = hichi.ScaleMapping(hichi.Axis.x, scaleCoef)


minInverseCoords = hichi.vector3d(-5, -5, -5)
maxInverseCoords = hichi.vector3d(5, 5, 5)

# transform of a point
point = hichi.vector3d(-1.0, -1.0, 0.0)
point2 = shiftMapping.getInverseCoords(rotationZMapping.getInverseCoords(point))
print("Direct mapping: point (%0.2f, %0.2f, %0.2f) -> (%0.2f, %0.2f, %0.2f)" % \
       (point.x, point.y, point.z, point2.x, point2.y, point2.z))
point3 = shiftMapping.getDirectCoords(rotationZMapping.getDirectCoords(point2))
print("Inverse mapping: point (%0.2f, %0.2f, %0.2f) -> (%0.2f, %0.2f, %0.2f)" % \
       (point2.x, point2.y, point2.z, point3.x, point3.y, point3.z))


gridSize = hichi.vector3d(64, 64, 64)

gridStep = hichi.vector3d(
            (maxInverseCoords.x - minInverseCoords.x) / gridSize.x,
            (maxInverseCoords.y - minInverseCoords.y) / gridSize.y,
            (maxInverseCoords.z - minInverseCoords.z) / gridSize.z
           )
           
timeStep = 0.1/hichi.c

# create computational grid
grid = hichi.PSATDGrid(gridSize, timeStep, minInverseCoords, gridStep)

# create shallow copy of computational grid
# gridMapping does direct transform when setters are called
# and inverse transform when getters are called
gridMapping1 = hichi.PSATDGridMapping(grid)
# set mapping to grid
gridMapping1.setMapping(rotationZMapping2)

# we can set more than one mapping
gridMapping2 = hichi.PSATDGridMapping(grid)
# set mapping to grid in correct order
gridMapping2.setMapping(scaleMapping)
gridMapping2.setMapping(rotationZMapping)
gridMapping2.setMapping(shiftMapping)

# create other grid with mappings which order is not correct
gridMapping3 = hichi.PSATDGridMapping(grid)
gridMapping3.setMapping(shiftMapping)
gridMapping3.setMapping(rotationZMapping)
gridMapping3.setMapping(scaleMapping)

# set fields to main grid
grid.setE(nullValue, fieldValue, nullValue)
grid.setB(nullValue, nullValue, fieldValue)

# create field solver
# periodic boundary conditions are default
fieldSolver = hichi.PSATDWithPoisson(grid)
# equivalent with
# fieldSolver = hichi.PSATDWithPoisson(gridMapping)

def updateData():
    fieldSolver.updateFields()


# --------- show -------------

import matplotlib.pyplot as plt
import matplotlib.animation as animation

N = 128
x = np.linspace(minInverseCoords.x, maxInverseCoords.x, N)
y = np.linspace(minInverseCoords.y, maxInverseCoords.y, N)


# get setted fields
def getFields(grid):
    global x, y, N
    Ey = np.zeros(shape=(N,N))
    for ix in range(N):
        for iy in range(N):
            coordXY = hichi.vector3d(x[ix], y[iy], 0.0)
            E = grid.getE(coordXY)
            Ey[N - iy - 1, ix] = E.y
    return Ey
    

Ey = getFields(grid)
transformedEy1 = getFields(gridMapping1)
transformedEy2 = getFields(gridMapping2)
transformedEy3 = getFields(gridMapping3)

fig, axes = plt.subplots(ncols=2, nrows=2)


# show setted field
im1 = axes[0, 0].imshow(Ey, cmap='RdBu', interpolation='none',
    extent=(minDirectCoords.x, maxDirectCoords.x, minDirectCoords.y, maxDirectCoords.y), animated = True)
fig.colorbar(im1, ax=axes[0, 0])
axes[0, 0].set_xlabel("x")
axes[0, 0].set_ylabel("y")
axes[0, 0].grid()
axes[0, 0].set_title("Setted Ey")


# show transformed field (rotation)
im2 = axes[0, 1].imshow(transformedEy1, cmap='RdBu', interpolation='none',
    extent=(minInverseCoords.x, maxInverseCoords.x, minInverseCoords.y, maxInverseCoords.y), animated = True)
fig.colorbar(im2, ax=axes[0, 1])
axes[0, 1].set_xlabel("x")
axes[0, 1].set_ylabel("y")
axes[0, 1].grid()
axes[0, 1].set_title("Transformed Ey \n rotation(pi/6)")


# show transformed field (rotation + shift)
im3 = axes[1, 0].imshow(transformedEy2, cmap='RdBu', interpolation='none',
    extent=(minInverseCoords.x, maxInverseCoords.x, minInverseCoords.y, maxInverseCoords.y), animated = True)
fig.colorbar(im3, ax=axes[1, 0])
axes[1, 0].set_xlabel("x")
axes[1, 0].set_ylabel("y")
axes[1, 0].grid()
axes[1, 0].set_title("Transformed Ey \n  scale(x, 1.5) + rotation(pi/3) \n + shift(-2, -2)")


# show transformed field (shift + rotation)
im4 = axes[1, 1].imshow(transformedEy3, cmap='RdBu', interpolation='none',
    extent=(minInverseCoords.x, maxInverseCoords.x, minInverseCoords.y, maxInverseCoords.y), animated = True)
fig.colorbar(im4, ax=axes[1, 1])
axes[1, 1].set_xlabel("x")
axes[1, 1].grid()
axes[1, 1].set_ylabel("y")
axes[1, 1].set_title("Transformed Ey \n shift(-2, -2) + rotation(pi/3) \n + scale(x, 1.5)")


def updatefig(*args):
    updateData()
    im1.set_array(getFields(grid))
    im2.set_array(getFields(gridMapping1))
    im3.set_array(getFields(gridMapping2))
    im4.set_array(getFields(gridMapping3))
    return im1, im2, im3, im4,
    
ani = animation.FuncAnimation(fig, updatefig, interval=10, blit=True)

fig.tight_layout()

plt.show()

