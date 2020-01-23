import sys
import os
import pyHiChi as hichi
import tight_focusing_fields as sphericalPulse
import tight_focusing_write_file as fileWriter
import math as ma
import numpy as np
import hichi_primitives

DIR_SCRIPT = os.path.abspath("./").replace("\\", "/")
DIR_RESULT = DIR_SCRIPT+"/results_theta/"

factor = 1
Nx = 32*factor
NxFull = 320*factor
Ny = factor*256
Nz = factor*256
gridSize = hichi.vector3d(Nx, Ny, Nz)

sphericalPulse.createSphericalPulseC()

timeStep = sphericalPulse.getDtCGS(16)
maxIter = 1

NFNumber = 20

D = 2*sphericalPulse.pulselength

mapping = hichi.TightFocusingMapping(sphericalPulse.R0, sphericalPulse.pulselength, D)

xMin = mapping.getxMin()
xMax = mapping.getxMax()

minCoords = hichi.vector3d(-20, -20, -20)
maxCoords = hichi.vector3d(20, 20, 20)

gridMinCoords = hichi.vector3d(xMin, minCoords.y, minCoords.z)
gridMaxCoords = hichi.vector3d(xMax, maxCoords.y, maxCoords.z)

gridStep = hichi.vector3d((gridMaxCoords.x - gridMinCoords.x) / gridSize.x, \
                          (gridMaxCoords.y - gridMinCoords.y) / gridSize.y, \
                          (gridMaxCoords.z - gridMinCoords.z) / gridSize.z)

thetaMax = ma.pi * 0.5 - 2*sphericalPulse.edgeSmoothingAngle
thetaMin = 2*sphericalPulse.edgeSmoothingAngle
thetaArr = np.arange(thetaMin, thetaMax, (thetaMax - thetaMin)/NFNumber)
F_number_arr = 1.0/(2.0 * np.tan(thetaArr))


#------- run -----------------------------


hichi_primitives.createDir(DIR_RESULT)  # DELETE OLD FILES IN DIR_RESULT!!!

ifCut = True
grid = hichi.PSATDGridMapping(gridSize, timeStep, gridMinCoords, gridStep, mapping, ifCut)

fieldSolver = hichi.PSATDWithPoisson(grid)

def updateFields():
    global mapping
    mapping.advanceTime(timeStep)
    fieldSolver.updateFields()

for theta in thetaArr:

    f_number = 1.0/(2.0 * ma.tan(theta))

    sphericalPulse.createSphericalPulseC(F_number_ = f_number)
    
    mapping.setTime(0.0)

    sphericalPulse.setField(grid)

    fileWriter.writeOX(grid, updateFields, minCoords, maxCoords, NxFull, maxIter=maxIter, dumpIter=maxIter,\
        fileName = "res_F_number_"+str(f_number)[:5]+"_iter_%d.csv", dirResult = DIR_RESULT, ifWriteZeroIter = False)
        
    
#-------- plot ------------------

import matplotlib.pyplot as plt
import matplotlib.animation as animation

arrPlotMax = []
arrPlotX = []

for f_number in F_number_arr:
 
    data = []
    with open(DIR_RESULT + "res_F_number_"+str(f_number)[:5]+"_iter_%d.csv" % maxIter, "r") as file:
        lines = file.readlines()
        for line in lines:
            data.append(float(line))
    
    m = max(data)
    arrPlotMax.append(m)
    arrPlotX.append(minCoords.x + data.index(m) * ((maxCoords.x - minCoords.x) / NxFull))

fig = plt.figure()
plt.plot(F_number_arr, arrPlotMax, "-o")
plt.xlabel('F_number')
plt.ylabel('max |E|')
plt.grid()
fig.tight_layout()
plt.savefig(DIR_RESULT + "max.png")

fig = plt.figure()
plt.plot(F_number_arr, arrPlotX, "-o")
plt.xlabel('F_number')
plt.ylabel('x_max')
plt.grid()
fig.tight_layout()
plt.savefig(DIR_RESULT + "xMax.png")

print(F_number_arr[arrPlotMax.index(max(arrPlotMax))], max(arrPlotMax), arrPlotX[arrPlotMax.index(max(arrPlotMax))])






