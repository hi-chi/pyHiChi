print('tightFocusing1')
import sys
sys.path.append("./../../build/src/pyHiChi/Release")
from math import *
import numpy as np
from pyHiChi import *
import pyHiChi as hichi
from hichi_primitives import *
import matplotlib.pyplot as plt
import tight_focusing_fields as sphericalPulse
import show

from tight_focusing_fields import R0
sizeX = (R0 + 5)*1e-4
sizeY = sizeX

ff = sphericalPulse.getFieldFuncs()
fieldDraft = hichi.analyticalField(ff[0][0], ff[0][1], ff[0][2], ff[1][0], ff[1][1], ff[1][2])
sphericalPulse = fieldDraft.init(0, 0.1*1e-4/lightVelocity)
#show.field(sphericalPulse, -sizeX, sizeX, -sizeY, sizeY, 256, 256, 30)
sizeZ = 0.7*sizeY
print('R0 = ', R0)
minCoords = hichi.vector3d(-(R0 + 3)*1e-4, -sizeZ, -sizeZ) # should be possible to write minCoords = 3*vector3d(-sizeX, -sizeY, -sizeZ) # we need to unify vector3d of primities and hichi
maxCoords = hichi.vector3d(-(R0 - 1)*1e-4, sizeZ, sizeZ)
#gridSize = hichi.vector3d(2*64, 2*64, 2*64)
factor = 1
gridSize = hichi.vector3d(factor*32, factor*256, factor*256)
fiedDraft = hichi.PSRTD(minCoords, maxCoords, gridSize, 1)
print('fieldDraft created.\n')
field1 = fiedDraft.init(0, 0.1*1e-4/lightVelocity)

MovingWindow = movingWindowCX(-(R0 + 3)*1e-4, -(R0 - 1)*1e-4, -(R0 + 1)*1e-4, field1)
field = MovingWindow(field1)
#field = field1

field.interpolateFrom(sphericalPulse)
dz = 2 * sizeZ / gridSize.z
def cond(x, y, z):
    return (dz / 2 - np.abs(z - dz / 2))
numberOfPictures = 64*8*3
import time
timeStart = time.time()
show.field(field, -sizeX, sizeX, -sizeY, sizeY, factor*256, factor*256, numberOfPictures, 8, func3ToC(cond))
#show.field(field, -sizeX, sizeX, -sizeY, sizeY, 512, 512, numberOfPictures, 8, True, func3ToC(cond))
totalTime = time.time() - timeStart
print('Total time = ', totalTime, '\n')
print("Time spent on field advance = %10.2f s, which is %10.2f %% of total time" % (field.getAdvanceTimeSpent(), 100*field.getAdvanceTimeSpent()/totalTime))
print("Performance = ", (1e+9)*(field.getAdvanceTimeSpent()/numberOfPictures)/(gridSize.x * gridSize.y * gridSize.z), " ns/cell.")
