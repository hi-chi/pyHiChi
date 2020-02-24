import sys
import os
sys.path.append("./../")
import pyHiChi as hichi
import tight_focusing_fields as sphericalPulse
import math as ma
import numpy as np
import hichi_primitives

DIR_SCRIPT = os.path.abspath("./").replace("\\", "/")
DIR_RESULT = DIR_SCRIPT+"/results_D/"

factor = 1.5
NxFull = int(320*factor)
Ny = int(256*factor)
Nz = int(256*factor)

minCoords = hichi.vector3d(-20, -20, -20)
maxCoords = hichi.vector3d(20, 20, 20)

gridStep = hichi.vector3d((maxCoords.x - minCoords.x) / NxFull, \
                          (maxCoords.y - minCoords.y) / Ny, \
                          (maxCoords.z - minCoords.z) / Nz)

sphericalPulse.createSphericalPulseC(F_number_ = 0.3)

timeStep = sphericalPulse.getDtCGS(sphericalPulse.R0)

Rmax = sphericalPulse.R0 + 0.5*sphericalPulse.pulselength
Rmin = sphericalPulse.R0 - 0.5*sphericalPulse.pulselength
angle = sphericalPulse.openingAngle + sphericalPulse.edgeSmoothingAngle
Dmin = -Rmin*np.cos(angle) + np.sqrt(Rmax**2 - (Rmin*np.sin(angle))**2)

NxMin = int(Dmin / gridStep.x) + 1
NxStep = 3 #2*factor
ND = 70
NxMax = NxMin + NxStep * ND

NxArr = np.arange(NxMin, NxMax, NxStep)
DArr = NxArr * gridStep.x


# ---------------------- run functions ----------------------


def updateFields(grid, mapping, fieldSolver):
    fieldSolver.setTimeStep(timeStep)
    mapping.advanceTime(timeStep)
    fieldSolver.updateFields()
    
    
def getFields(grid):
    x = np.arange(minCoords.x, maxCoords.x, (maxCoords.x - minCoords.x)/NxFull)
    y = np.arange(minCoords.y, maxCoords.y, (maxCoords.y - minCoords.y)/Ny)
    field = np.zeros(shape=(NxFull,Ny))
    for ix in range(NxFull):
        for iy in range(Ny):
            coordXY = hichi.vector3d(x[ix], y[iy], 0.0)
            E = grid.getE(coordXY)
            field[ix, iy] = E.norm()
    return field


def run(Nx):

    D = Nx * gridStep.x

    mapping = hichi.TightFocusingMapping(sphericalPulse.R0, sphericalPulse.pulselength, D)
    #mapping.setIfCut(False)

    xMin = mapping.getxMin()
    xMax = mapping.getxMax()
    
    gridMinCoords = hichi.vector3d(xMin, minCoords.y, minCoords.z)
    gridMaxCoords = hichi.vector3d(xMax, maxCoords.y, maxCoords.z)
    
    gridSize = hichi.vector3d(Nx, Ny, Nz)
    
    grid = hichi.PSATDGridMapping(gridSize, timeStep, gridMinCoords, gridStep)
    grid.setMapping(mapping)
    
    fieldSolver = hichi.PSATDWithPoisson(grid)
    
    mapping.setTime(0.0)

    sphericalPulse.setField(grid)
    
    updateFields(grid, mapping, fieldSolver)
    
    return getFields(grid)
    
    
def computeError(resBand, resFull):
    m = 0.0
    for iy in range(Ny):
        for ix in range(NxFull):
            diff = abs(resBand[ix, iy] - resFull[ix, iy])
            m = diff if m < diff else m
    return m


# ---------------------- run -----------------------


hichi_primitives.createDir(DIR_RESULT)

resDiff = []

resFull = run(NxFull)


for Nx in NxArr:

    print("\r%d" % Nx, end="")

    resBand = run(Nx)
    
    resDiff.append(computeError(resBand, resFull))
    

# ---------------------- save -----------------------


with open(DIR_RESULT + "/res_D.csv", "w") as file:
    for i in range(len(resDiff)):
        file.write("%d;%0.15f;%0.15f\n" % (NxArr[i], DArr[i], resDiff[i]))
    
    
# ---------------------- read -----------------------

resDiff = []
with open(DIR_RESULT + "/res_D.csv", "r") as file:
    for line in file.readlines():
        arr = [float(elem) for elem in line.split(";")]
        resDiff.append(arr[2])
        

# ---------------------- plot -----------------------


import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation

matplotlib.rcParams.update({"font.size" : 17})

fig = plt.figure()
plt.plot(DArr, resDiff, "-3b")
plt.xlabel('D')
plt.ylabel('error')
plt.grid()
fig.tight_layout()
plt.savefig(DIR_RESULT + "/error_D.png")

