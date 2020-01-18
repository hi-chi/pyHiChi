import sys
import os
sys.path.append("../../build/src/pyHiChi/Release")
import pyHiChi as hichi
import tight_focusing_fields as sphericalPulse
import tight_focusing_write_file as fileWriter
import math as ma
import numpy as np
import hichi_primitives

DIR_SCRIPT = os.path.abspath("./").replace("\\", "/")
DIR_RESULT = DIR_SCRIPT+"/results/"

factor = 1
Nx = 32*factor
NxFull = 320*factor
gridSize = hichi.vector3d(Nx, factor*256, factor*256)

timeStep = 1*sphericalPulse.wavelength/hichi.c
maxIter = 16

NFNumber = 20

D = 2*sphericalPulse.pulseLength

mapping = hichi.TightFocusingMapping(sphericalPulse.R0, sphericalPulse.pulseLength, D)

xMin = mapping.getxMin()
xMax = mapping.getxMax()

minCoords = hichi.vector3d(-20e-4, -20e-4, -20e-4)
maxCoords = hichi.vector3d(20e-4, 20e-4, 20e-4)

gridMinCoords = hichi.vector3d(xMin, minCoords.y, minCoords.z)
gridMaxCoords = hichi.vector3d(xMax, maxCoords.y, maxCoords.z)

gridStep = hichi.vector3d((gridMaxCoords.x - gridMinCoords.x) / gridSize.x, \
                          (gridMaxCoords.y - gridMinCoords.y) / gridSize.y, \
                          (gridMaxCoords.z - gridMinCoords.z) / gridSize.z)

thetaMax = ma.pi * 0.5 - 2*sphericalPulse.edgeSmoothingAngle
thetaMin = 2*sphericalPulse.edgeSmoothingAngle
thetaArr = np.arange(thetaMin, thetaMax, (thetaMax - thetaMin)/NFNumber)
F_number_min = 1.0/(2.0 * ma.tan(thetaMax))
F_number_max = 1.0/(2.0 * ma.tan(thetaMin))
F_number_arr = 1.0/(2.0 * np.tan(thetaArr))
#F_number_arr = np.arange(F_number_min, F_number_max, (F_number_max - F_number_min)/NFNumber)


#------- run -----------------------------


hichi_primitives.createDir(DIR_RESULT)

ifCut = True
grid = hichi.PSATDGridMapping(gridSize, timeStep, gridMinCoords, gridStep, mapping, ifCut)

fieldSolver = hichi.PSATDWithPoisson(grid)

def updateFields():
    global mapping
    mapping.advanceTime(timeStep)
    fieldSolver.updateFields()

for f_number in F_number_arr:

    sphericalPulse.createSphericalPulseC(F_number_ = f_number)
    
    mapping.setTime(0.0)

    sphericalPulse.setField(grid)
    
    fileWriter.writeOX(grid, updateFields, minCoords, maxCoords, NxFull, maxIter=maxIter, dumpIter=maxIter,\
        fileName = "res_F_number_"+str(f_number)+"_iter_%d.csv", dirResult = DIR_RESULT, ifWriteZeroIter = False)
    
    
    
#-------- plot ------------------

# import matplotlib.pyplot as plt
# import matplotlib.animation as animation

# arrPlotMax = []
# arrPlotX = []

# for f_number in F_number_arr:
 
    # data = []
    # with open(DIR_RESULT + "res_F_number_"+str(f_number)+"_iter_%d.csv" % maxIter, "r") as file:
        # lines = file.readlines()
        # for line in lines:
            # data.append(float(line))
    
    # m = max(data)
    # arrPlotMax.append(m)
    # arrPlotX.append(minCoords.x + data.index(m) * ((maxCoords.x - minCoords.x) / NxFull))

# fig = plt.figure()
# plt.plot(F_number_arr, arrPlotMax, "*-")
# plt.savefig(DIR_RESULT + "max.png")

# fig = plt.figure()
# plt.plot(F_number_arr, arrPlotX, "*-")
# plt.savefig(DIR_RESULT + "xMax.png")





