import sys
sys.path.append("../../build/src/pyHiChi/Release")
import pyHiChi as hichi
import tight_focusing_fields as sphericalPulse
import tight_focusing_write_file as fileWriter
import math as ma
import numpy as np
import hichi_primitives

# import matplotlib.pyplot as plt
# import matplotlib.animation as animation

hichi_primitives.createDir("./results/")

factor = 1
Nx = 32
gridSize = hichi.vector3d(factor*Nx, factor*256, factor*256)

timeStep = 0.2*sphericalPulse.wavelength/hichi.c
maxIter = 81

NFNumber = 20

D = 2*sphericalPulse.pulseLength

mapping = hichi.TightFocusingMapping(sphericalPulse.R0, sphericalPulse.pulseLength, D)

xMin = mapping.getxMin()
xMax = mapping.getxMax()

widthCoeff = 2
coordMax = widthCoeff*sphericalPulse.R0
minCoords = hichi.vector3d(-coordMax, -coordMax, -coordMax)
maxCoords = hichi.vector3d(coordMax, coordMax, coordMax)

gridMinCoords = hichi.vector3d(xMin, minCoords.y, minCoords.z)
gridMaxCoords = hichi.vector3d(xMax, maxCoords.y, maxCoords.z)

gridStep = hichi.vector3d((gridMaxCoords.x - gridMinCoords.x) / gridSize.x, \
                          (gridMaxCoords.y - gridMinCoords.y) / gridSize.y, \
                          (gridMaxCoords.z - gridMinCoords.z) / gridSize.z)

ifCut = True
grid = hichi.PSATDGridMapping(gridSize, timeStep, gridMinCoords, gridStep, mapping, ifCut)

fieldSolver = hichi.PSATDWithPoisson(grid)

def updateFields():
    global mapping
    mapping.advanceTime(timeStep)
    fieldSolver.updateFields()

thetaMax = ma.pi * 0.5 - sphericalPulse.edgeSmoothingAngle
thetaMin = 2*sphericalPulse.edgeSmoothingAngle
F_number_min = 1.0/(2.0 * ma.tan(thetaMax))
F_number_max = 1.0/(2.0 * ma.tan(thetaMin))
F_number_arr = np.arange(F_number_min, F_number_max, (F_number_max - F_number_min)/NFNumber)

# arrPlotMax = []
# arrPlotX = []

for f_number in F_number_arr:

    sphericalPulse.createSphericalPulse(F_number_ = f_number)
    
    mapping.setTime(0.0)

    fieldFuncs = sphericalPulse.getFieldFuncs()
    grid.setE(fieldFuncs[0][0], fieldFuncs[0][1], fieldFuncs[0][2])
    grid.setB(fieldFuncs[1][0], fieldFuncs[1][1], fieldFuncs[1][2])
    
    fileWriter.writeOX(grid, updateFields, minCoords, maxCoords, Nx, maxIter=maxIter, dumpIter=maxIter,\
        fileName = "res_F_number_"+str(f_number)+"_iter_%d.csv", dirResult = "./results/", ifWriteZeroIter = False)
    
    # data = []
    # with open("./results/res_F_number_"+str(f_number)+"_iter_%d.csv" % maxIter, "r") as file:
        # lines = file.readlines()
        # for line in lines:
            # data.append(float(line))
    
    #m = max(data)
    #arrPlotMax.append(max(data))
    #arrPlotX.append(minCoords.x + data.index(m) * (maxCoords.x - minCoords.y) / Nx)

# fig = plt.figure()
# plt.plot(F_number_arr, arrPlotMax, "*-")
# plt.savefig("max.png")

# fig = plt.figure()
# plt.plot(F_number_arr, arrPlotX, "*-")
# plt.savefig("xMax.png")


