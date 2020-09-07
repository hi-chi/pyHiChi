import sys
sys.path.append("../bin")
import pyHiChi as hichi
import numpy as np


minCoords = hichi.vector3d(-10, -10, 0.0)
maxCoords = hichi.vector3d(10, 10, 0.0)

    
def nullValue(x, y, z):
    return 0.0


# ------- create the first pulse ---------------

def fieldValue1(x, y, z):
    return np.exp(-x**2-y**2)*np.sin(3*x)  # omega=5
    
gridSize1 = hichi.vector3d(128, 128, 1)
d1 = (maxCoords.x - minCoords.x) / gridSize1.x
gridStep1 = hichi.vector3d(d1, d1, d1)
timeStep1 = 0.05/hichi.c

field1 = hichi.PSATDField(gridSize1, minCoords, gridStep1, timeStep1)
field1.setE(nullValue, fieldValue1, nullValue)
field1.setB(nullValue, nullValue, fieldValue1)
field1.convertFieldsPoissonEquation()


# ------- create the second pulse --------------

def fieldValue2(x, y, z):
    return np.exp(-x**2-y**2)*np.sin(6*x)  # omega=10
    
gridSize2 = hichi.vector3d(256, 256, 1)
d2 = (maxCoords.x - minCoords.x) / gridSize2.x
gridStep2 = hichi.vector3d(d2, d2, d2)
timeStep2 = 0.025/hichi.c  # PSTD Courant condition is taken into account

field2 = hichi.PSTDField(gridSize2, minCoords, gridStep2, timeStep2)
field2.setE(nullValue, fieldValue2, nullValue)
field2.setB(nullValue, nullValue, fieldValue2)


# ------ transform fields ---------------------

rotationZMapping = hichi.RotationMapping(hichi.Axis.z, hichi.pi/3)
shiftMapping1 = hichi.ShiftMapping(hichi.vector3d(-6.5, -4.0, 0.0))
shiftMapping2 = hichi.ShiftMapping(hichi.vector3d(5.0, -4.0, 0.0))
shiftMapping3 = hichi.ShiftMapping(hichi.vector3d(-7.0, 7.0, 0.0))

resField = (0.5*field2.applyMapping(shiftMapping1) + field1 + (field2*0.8).applyMapping(shiftMapping2))\
                .applyMapping(rotationZMapping) + field1.applyMapping(shiftMapping3)


# ------ update fields ------------------------

def updateData():  # dt = 0.05/hichi.c
    field1.updateFields()  # dt = 0.05/hichi.c
    field2.updateFields()  # dt = 0.025/hichi.c
    field2.updateFields()  # dt = 0.025/hichi.c
    
    
# --------- show -------------

import matplotlib.pyplot as plt
import matplotlib.animation as animation

N = 128
x = np.linspace(minCoords.x, maxCoords.x, N)
y = np.linspace(minCoords.y, maxCoords.y, N)

def getFields(field):
    global x, y, N
    res = np.zeros(shape=(N,N))
    for ix in range(N):
        for iy in range(N):
            coordXY = hichi.vector3d(x[ix], y[iy], 0.0)
            res[N - iy - 1, ix] = field.getE(coordXY).norm()
    return res

fig = plt.figure()
ax = fig.add_subplot(1,1,1)

im = ax.imshow(getFields(resField), cmap='RdBu', interpolation='none',
    extent=(minCoords.x, maxCoords.x, minCoords.y, maxCoords.y), animated = True)
fig.colorbar(im, ax=ax)
ax.set_xlabel("x")
ax.set_ylabel("y")

def updatefig(*args):
    updateData()
    im.set_array(getFields(resField))
    return im,
    
ani = animation.FuncAnimation(fig, updatefig, interval=10, blit=True)

fig.tight_layout()

plt.show()


