import sys
import os
import pyHiChi as hichi
import tight_focusing_fields as sphericalPulse
import tight_focusing_write_file as fileWriter
import math as ma
import numpy as np
import hichi_primitives

LAUNCHE_FOR_SMALL_GRID = "small"
LAUNCHE_FOR_LARGE_GRID = "large"

#LAUNCHE = LAUNCHE_FOR_SMALL_GRID
LAUNCHE = LAUNCHE_FOR_LARGE_GRID

DIR_SCRIPT = os.path.abspath("./").replace("\\", "/")

DIR_RESULT_SMALL = DIR_SCRIPT+"/results_theta_%s/" % LAUNCHE_FOR_SMALL_GRID
DIR_RESULT_LARGE = DIR_SCRIPT+"/results_theta_%s/" % LAUNCHE_FOR_LARGE_GRID

DIR_RESULT = DIR_RESULT_SMALL if LAUNCHE == LAUNCHE_FOR_SMALL_GRID else DIR_RESULT_LARGE

factor = 1 if LAUNCHE == LAUNCHE_FOR_SMALL_GRID else 6
Nx = int(32*factor)
NxFull = int(320*factor)
Ny = int(factor*256)
Nz = int(factor*256)
gridSize = hichi.vector3d(Nx, Ny, Nz)

sphericalPulse.createSphericalPulseC()

NFNumber = 20

D = 2*sphericalPulse.pulselength

mapping = hichi.TightFocusingMapping(sphericalPulse.R0, sphericalPulse.pulselength, D)

xMin = mapping.getxMin()
xMax = mapping.getxMax()

minCoords = hichi.vector3d(-20, -20, -20)
maxCoords = hichi.vector3d(20, 20, 20)
dx = (maxCoords.x-minCoords.x)/NxFull

gridMinCoords = hichi.vector3d(xMin, minCoords.y, minCoords.z)
gridMaxCoords = hichi.vector3d(xMax, maxCoords.y, maxCoords.z)

gridStep = hichi.vector3d((gridMaxCoords.x - gridMinCoords.x) / gridSize.x, \
                          (gridMaxCoords.y - gridMinCoords.y) / gridSize.y, \
                          (gridMaxCoords.z - gridMinCoords.z) / gridSize.z)

thetaMax = ma.pi * 0.5 - 2*sphericalPulse.edgeSmoothingAngle
thetaMin = 2*sphericalPulse.edgeSmoothingAngle
#thetaArr = np.arange(thetaMin, thetaMax, (thetaMax - thetaMin)/NFNumber)
def FNumber(theta):
    return 1.0/(2.0 * np.tan(theta))
F_number_arr = np.arange(FNumber(thetaMin), FNumber(thetaMax), (FNumber(thetaMax) - FNumber(thetaMin))/NFNumber)


#------- run -------------------------------


hichi_primitives.createDir(DIR_RESULT)  # DELETE OLD FILES IN DIR_RESULT!!!

grid = hichi.PSATDGridMapping(gridSize, 0.0, gridMinCoords, gridStep)
grid.setMapping(mapping)

fieldSolver = hichi.PSATDWithPoisson(grid)

def updateFields(timeStep):
    fieldSolver.setTimeStep(sphericalPulse.getDtCGS(timeStep))
    mapping.advanceTime(sphericalPulse.getDtCGS(timeStep))
    fieldSolver.updateFields()
      
def getFields():
    x = np.arange(minCoords.x, maxCoords.x, dx)
    field = np.zeros(shape=(NxFull))
    for ix in range(NxFull):
        coordXY = hichi.vector3d(x[ix], 0.0, 0.0)
        E = grid.getE(coordXY)
        field[ix] = E.norm()
    return field

  
# run for small grid
def runSmall(f_number):
    timeStep = sphericalPulse.R0
    updateFields(timeStep)
    
    timeStep = -0.1
    maxIter = 180
    
    with open(DIR_RESULT + "/res_%0.3f.csv" % f_number, "w") as file:
        for iter in range(maxIter+1):
            field = getFields()
            m = field.max()
            xmax = minCoords.x + dx * np.where(field == m)[0]
            file.write("%0.15f;%0.15f\n" %(m, xmax))
            updateFields(timeStep)


# run for large grid
def runLarge(f_number, xMiddle):
    eps = 1.5
    xmin = max(xMiddle - eps, -sphericalPulse.R0)
    xmax = xMiddle + eps
    tstart = xmin + sphericalPulse.R0
    tend   = xmax + sphericalPulse.R0
    
    timeStep = tstart
    updateFields(timeStep)
    
    timeStep = 0.02
    maxIter = int((tend - tstart)/timeStep) + 1
    print(xMiddle, xmin, xmax, tstart, tend, timeStep, maxIter)
    
    with open(DIR_RESULT + "/res_%0.3f.csv" % f_number, "w") as file:
        for iter in range(maxIter+1):
            field = getFields()
            m = field.max()
            xmax = minCoords.x + dx * np.where(field == m)[0]
            file.write("%0.15f;%0.15f\n" % (m, xmax))
            updateFields(timeStep)


def readSmallResults():
    res = {}
    with open(DIR_RESULT_SMALL + "/res.csv", 'r') as file:
        for line in file.readlines():
            l = line.split(";")
            res[l[0]] = (float(l[1]), float(l[2]))
    return res


if (LAUNCHE == LAUNCHE_FOR_LARGE_GRID):
    smallResults = readSmallResults()


for f_number in F_number_arr:

    sphericalPulse.createSphericalPulseC(F_number_ = f_number)
    
    mapping.setTime(0.0)

    sphericalPulse.setField(grid)
    
    if (LAUNCHE == LAUNCHE_FOR_SMALL_GRID):
        runSmall(f_number)
    else:
        runLarge(f_number, smallResults["%0.15f" % f_number][1])
        
    
#-------- plot ------------------

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation

matplotlib.rcParams.update({"font.size" : 17})

arrPlotMax = []
arrPlotX = []

for f_number in F_number_arr:
 
    data = [[],[]]
    with open(DIR_RESULT + ("res_%0.3f" % f_number)+".csv", "r") as file:
        lines = file.readlines()
        for line in lines:
            arr = [float(elem) for elem in line.split(";")]
            data[0].append(arr[0])
            data[1].append(arr[1])
    
    m = max(data[0])
    ind = data[0].index(m)
    arrPlotMax.append(m)
    arrPlotX.append(data[1][ind])
    
with open(DIR_RESULT + "/res.csv", "w") as file:
    for f_number, max, xmax in zip(F_number_arr, arrPlotMax, arrPlotX):
        file.write("%0.15f;%0.15f;%0.15f\n" % (f_number, max, xmax))

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


# import matplotlib
# import matplotlib.pyplot as plt
# import matplotlib.animation as animation

# matplotlib.rcParams.update({"font.size" : 17})


# for f_number in F_number_arr:
 
    # data = [[],[]]
    # with open(DIR_RESULT + ("res_%0.3f" % f_number)+".csv", "r") as file:
        # lines = file.readlines()
        # for line in lines:
            # arr = [float(elem) for elem in line.split(";")]
            # data[0].append(arr[0])
            # data[1].append(arr[1])
    
    # fig = plt.figure()
    # plt.plot(data[1], data[0], "-*")
    # plt.savefig(DIR_RESULT + "%0.3f.png" % f_number)


