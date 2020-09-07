import sys
sys.path.append("../bin")
import pyHiChi as hichi
import numpy as np


def fieldValue(x, y, z):
    return np.exp(-x**2-y**2-z**2)*np.sin(5*x)
    
def nullValue(x, y, z):
    return 0.0
    
minCoords = hichi.vector3d(-5, -5, -5)
maxCoords = hichi.vector3d(5, 5, 5)

gridSize = hichi.vector3d(64, 64, 64)
gridStep = (maxCoords - minCoords) / gridSize  
timeStep = 0.1/hichi.c

field = hichi.PSATDField(gridSize, minCoords, gridStep, timeStep)
field.setE(nullValue, fieldValue, nullValue)
field.setB(nullValue, nullValue, fieldValue)
field.convertFieldsPoissonEquation()


# create the first transformed field
# rotation
def getRotatedField(field):
    mapping = hichi.RotationMapping(hichi.Axis.z, -np.pi/6)  # local object
    return field.applyMapping(mapping)  # mapping is inaccessible but not destroyed

field1 = getRotatedField(field)  


# create some other mappings  
angle = np.pi/3
rotationZMapping = hichi.RotationMapping(hichi.Axis.z, angle)

shift = hichi.vector3d(-2.0, -2.0, 0.0)
shiftMapping = hichi.ShiftMapping(shift)

scaleCoef = 1.5
scaleMapping = hichi.ScaleMapping(hichi.Axis.x, scaleCoef)


# transform a point (rotation -> shift)
point = hichi.vector3d(-1.0, -1.0, 0.0)
point2 = shiftMapping.getInverseCoords(rotationZMapping.getInverseCoords(point))
print("Direct mapping: point (%0.2f, %0.2f, %0.2f) -> (%0.2f, %0.2f, %0.2f)" % \
       (point.x, point.y, point.z, point2.x, point2.y, point2.z))
point3 = shiftMapping.getDirectCoords(rotationZMapping.getDirectCoords(point2))
print("Inverse mapping: point (%0.2f, %0.2f, %0.2f) -> (%0.2f, %0.2f, %0.2f)" % \
       (point2.x, point2.y, point2.z, point3.x, point3.y, point3.z))


# create the second transformed field
# scale -> rotation -> shift
field2Tmp1 = field.applyMapping(scaleMapping)
field2Tmp2 = field2Tmp1.applyMapping(rotationZMapping)
field2 = field2Tmp2.applyMapping(shiftMapping)


# create the third transformed field
# shift -> rotation -> scale
field3 = field.applyMapping(shiftMapping)\
              .applyMapping(rotationZMapping)\
              .applyMapping(scaleMapping)
# intermediate fields are inaccesible but not destroyed

def updateData():
    field.updateFields()
    # == field1.updateFields()
    # == field2.updateFields()
    # == field3.updateFields()
    
    
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

fig, axes = plt.subplots(ncols=2, nrows=2)

# show setted field
im1 = axes[0, 0].imshow(getFields(field), cmap='RdBu', interpolation='none',
    extent=(minCoords.x, maxCoords.x, minCoords.y, maxCoords.y), animated = True)
fig.colorbar(im1, ax=axes[0, 0])
axes[0, 0].set_xlabel("x")
axes[0, 0].set_ylabel("y")
axes[0, 0].set_title("Original field")


# show transformed field (rotation)
im2 = axes[0, 1].imshow(getFields(field1), cmap='RdBu', interpolation='none',
    extent=(minCoords.x, maxCoords.x, minCoords.y, maxCoords.y), animated = True)
fig.colorbar(im2, ax=axes[0, 1])
axes[0, 1].set_xlabel("x")
axes[0, 1].set_ylabel("y")
axes[0, 1].set_title("Transformed field \n rotation")


# show transformed field (scale -> rotation -> shift)
im3 = axes[1, 0].imshow(getFields(field2), cmap='RdBu', interpolation='none',
    extent=(minCoords.x, maxCoords.x, minCoords.y, maxCoords.y), animated = True)
fig.colorbar(im3, ax=axes[1, 0])
axes[1, 0].set_xlabel("x")
axes[1, 0].set_ylabel("y")
axes[1, 0].set_title("Transformed field \n scale -> rotation -> shift")


# show transformed field (shift + rotation)
im4 = axes[1, 1].imshow(getFields(field3), cmap='RdBu', interpolation='none',
    extent=(minCoords.x, maxCoords.x, minCoords.y, maxCoords.y), animated = True)
fig.colorbar(im4, ax=axes[1, 1])
axes[1, 1].set_xlabel("x")
axes[1, 1].set_ylabel("y")
axes[1, 1].set_title("Transformed field \n shift ->rotation -> scale")


def updatefig(*args):
    updateData()
    im1.set_array(getFields(field))
    im2.set_array(getFields(field1))
    im3.set_array(getFields(field2))
    im4.set_array(getFields(field3))
    return im1, im2, im3, im4,
    
ani = animation.FuncAnimation(fig, updatefig, interval=10, blit=True)

fig.tight_layout()

plt.show()
