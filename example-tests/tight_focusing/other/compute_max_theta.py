import sys
import os
sys.path.append("./../")
import pyHiChi as hichi
import tight_focusing_fields as sphericalPulse
#import tight_focusing_write_file as fileWriter
import math as ma
import numpy as np
import hichi_primitives

DIR_SCRIPT = os.path.abspath("./").replace("\\", "/")
DIR_RESULT = DIR_SCRIPT+"/results_theta/"

factor = 1.0
Nx = int(32*factor)
NxFull = int(320*factor)
Ny = int(factor*256)
Nz = int(factor*256)
gridSize = hichi.vector3d(Nx, Ny, Nz)

sphericalPulse.createSphericalPulsePython()

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

FNumber = lambda theta: 1.0/(2.0 * np.tan(theta))                          
thetaMax = ma.pi * 0.5 - 2*sphericalPulse.edgeSmoothingAngle
thetaMin = 2*sphericalPulse.edgeSmoothingAngle
NFNumber = 15
F_number_arr = np.arange(FNumber(thetaMin), FNumber(thetaMax), (FNumber(thetaMax) - FNumber(thetaMin))/NFNumber)


# ---------------------- run functions ----------------------


def updateFields(timeStep):
    fieldSolver.setTimeStep(sphericalPulse.getDtCGS(timeStep))
    mapping.advanceTime(sphericalPulse.getDtCGS(timeStep))
    fieldSolver.updateFields()


def getFields():
    x = np.arange(minCoords.x, maxCoords.x, dx)
    field = np.zeros(shape=(NxFull))
    coord = np.zeros(shape=(NxFull))
    for ix in range(NxFull):
        coordXY = hichi.vector3d(x[ix], 0.0, 0.0)
        E = grid.getE(coordXY)
        field[ix] = E.norm()
        coord[ix] = x[ix]
    return field, coord


def runInterval(timeStart, timeEnd, nIter):
    resMax = {}
    timeStep = (timeEnd - timeStart)/nIter
    time = timeStart
    for iter in range(nIter):
        field, coord = getFields()
        m = field.max()
        xm = coord[np.where(field == m)[0]]
        resMax[time] = (xm, m)
        time += timeStep
        updateFields(timeStep)
    return resMax


# nIter = 202, dt = 0.01
def run():
    res = []
    
    # dt = 1
    timeStart = 0
    timeEnd = 16
    nIter = 32
    resTmp = runInterval(timeStart, timeEnd, nIter)
    m = max(resTmp.items(), key = lambda m: m[1][1])
    timeMax = m[0]
    res.extend(resTmp.items())
    res.sort()
    
    # dt = 0.1
    timeStart = timeMax - 1.0
    updateFields(timeStart - timeEnd)
    timeEnd = timeMax + 1.0
    nIter = 20
    resTmp = runInterval(timeStart, timeEnd, nIter)
    m = max(resTmp.items(), key = lambda m: m[1][1])
    timeMax = m[0]
    res.extend(resTmp.items())
    res.sort()
    
    # dt = 0.01
    timeStart = timeMax - 0.75
    updateFields(timeStart - timeEnd)
    timeEnd = timeMax + 0.75
    nIter = 50
    resTmp = runInterval(timeStart, timeEnd, nIter)
    m = max(resTmp.items(), key = lambda m: m[1][1])
    timeMax = m[0]
    res.extend(resTmp.items())
    res.sort()
    
    return res


# nIter = 80, dt = 0.01
def runUsingPreviousResults(fileName, f_number):
    res = []
    timeMax = None
    timeStart = 0.0
    timeEnd = 0.0
    
    with open(fileName, "r") as file:
        for line in file.readlines():
            arr = line.split(";")
            if (arr[0] == "%0.15f" % f_number):
                timeMax = float(arr[1])
    
    # dt = 0.01
    timeStart = timeMax - 0.4
    updateFields(timeStart - timeEnd)
    timeEnd = timeMax + 0.4
    nIter = 80
    resTmp = runInterval(timeStart, timeEnd, nIter)
    m = max(resTmp.items(), key = lambda m: m[1][1])
    timeMax = m[0]
    res.extend(resTmp.items())
    res.sort()
    
    return res
    

# nIter = 160, dt = 0.1
def runConst():
    res = []
    
    # dt = 0.1
    timeStart = 0
    timeEnd = 16
    nIter = 160
    resTmp = runInterval(timeStart, timeEnd, nIter)
    m = max(resTmp.items(), key = lambda m: m[1][1])
    timeMax = m[0]
    res.extend(resTmp.items())
    res.sort()
    
    return res


# ---------------------- run ----------------------


hichi_primitives.createDir(DIR_RESULT)  # DELETE OLD FILES IN DIR_RESULT!!!

grid = hichi.PSATDGridMapping(gridSize, 0.0, gridMinCoords, gridStep)
grid.setMapping(mapping)

fieldSolver = hichi.PSATD(grid)


for f_number in F_number_arr:

    sphericalPulse.createSphericalPulsePython(F_number_ = f_number)
    
    mapping.setTime(0.0)

    sphericalPulse.setField(grid)
    fieldSolver.convertFieldsPoissonEquation()
    
    res = run()
    
    with open(DIR_RESULT + "/f_number_%0.3f.csv" % f_number, "w") as file:
        for time, value in res:
            file.write("%0.15f;%0.15f;%0.15f\n" % (time, value[0], value[1])) 


# for f_number in F_number_arr:

    # sphericalPulse.createSphericalPulseC(F_number_ = f_number)
    
    # mapping.setTime(0.0)

    # sphericalPulse.setField(grid)
    # fieldSolver.convertFieldsPoissonEquation()
    
    # res = runUsingPreviousResults(DIR_SCRIPT + "/res.csv", f_number)
    
    # with open(DIR_RESULT + "/f_number_%0.3f.csv" % f_number, "w") as file:
        # for time, value in res:
            # file.write("%0.15f;%0.15f;%0.15f\n" % (time, value[0], value[1]))


# # for f_number in F_number_arr:

    # # sphericalPulse.createSphericalPulseC(F_number_ = f_number)
    
    # # mapping.setTime(0.0)

    # # sphericalPulse.setField(grid)
    # # fieldSolver.convertFieldsPoissonEquation()
    
    # # res = runConst()
    
    # # with open(DIR_RESULT + "/f_number_const_%0.3f.csv" % f_number, "w") as file:
        # # for time, value in res:
            # # file.write("%0.15f;%0.15f;%0.15f\n" % (time, value[0], value[1])) 


# ---------------------- plot ----------------------



import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation

matplotlib.rcParams.update({"font.size" : 17})

arrPlotMax = []
arrPlotX = []
arrPlotTime = []

for f_number in F_number_arr:
 
    data = {"time":[], "xmax":[], "max":[]}
    with open(DIR_RESULT + "/f_number_%0.3f.csv" % f_number, "r") as file:
        lines = file.readlines()
        for line in lines:
            arr = [float(elem) for elem in line.split(";")]
            data["time"].append(arr[0])
            data["xmax"].append(arr[1])
            data["max"].append(arr[2])
    
    m = max(data["max"])
    ind = data["max"].index(m)
    arrPlotMax.append(m)
    arrPlotX.append(data["xmax"][ind])
    arrPlotTime.append(data["time"][ind])

    
with open(DIR_RESULT + "/res.csv", "w") as file:
    for F_number, time, xm, m in zip(F_number_arr, arrPlotTime, arrPlotX, arrPlotMax):
        file.write("%0.15f;%0.15f;%0.15f;%0.15f\n" % (F_number, time, xm, m))


# arrPlotTimeConst = []
# arrPlotMaxConst = []
# arrPlotXConst = []

# for f_number in F_number_arr:
 
    # data = {"time":[], "xmax":[], "max":[]}
    # with open(DIR_RESULT + "/f_number_const_%0.3f.csv" % f_number, "r") as file:
        # lines = file.readlines()
        # for line in lines:
            # arr = [float(elem) for elem in line.split(";")]
            # data["time"].append(arr[0])
            # data["xmax"].append(arr[1])
            # data["max"].append(arr[2])
    
    # m = max(data["max"])
    # ind = data["max"].index(m)
    # arrPlotMaxConst.append(m)
    # arrPlotXConst.append(data["xmax"][ind])
    # arrPlotTimeConst.append(data["time"][ind])

    
# with open(DIR_RESULT + "/res_const.csv", "w") as file:
    # for F_number, time, xm, m in zip(F_number_arr, arrPlotTimeConst, arrPlotXConst, arrPlotMaxConst):
        # file.write("%0.15f;%0.15f;%0.15f;%0.15f\n" % (F_number, time, xm, m))


fig = plt.figure()
#plt.plot(F_number_arr, arrPlotMaxConst, "-ob")
plt.plot(F_number_arr, arrPlotMax, "-ob")
plt.xlabel('F_number')
plt.ylabel('max |E|')
plt.grid()
fig.tight_layout()
plt.savefig(DIR_RESULT + "/max.png")

fig = plt.figure()
#plt.plot(F_number_arr, arrPlotXConst, "-ob")
plt.plot(F_number_arr, arrPlotX, "-ob")
plt.xlabel('F_number')
plt.ylabel('x_max')
plt.grid()
fig.tight_layout()
plt.savefig(DIR_RESULT + "/xMax.png")


for f_number in F_number_arr:

    fig = plt.figure()
 
    # data = {"time":[], "xmax":[], "max":[]}
    # with open(DIR_RESULT + "/f_number_const_%0.3f.csv" % f_number, "r") as file:
        # lines = file.readlines()
        # for line in lines:
            # arr = [float(elem) for elem in line.split(";")]
            # data["time"].append(arr[0])
            # data["xmax"].append(arr[1])
            # data["max"].append(arr[2])
    
    # plt.plot(data["time"], data["max"], "-b")
    
    
    data = {"time":[], "xmax":[], "max":[]}
    with open(DIR_RESULT + "/f_number_%0.3f.csv" % f_number, "r") as file:
        lines = file.readlines()
        for line in lines:
            arr = [float(elem) for elem in line.split(";")]
            data["time"].append(arr[0])
            data["xmax"].append(arr[1])
            data["max"].append(arr[2])
    
    plt.plot(data["time"], data["max"], "-r")
    
    plt.savefig(DIR_RESULT + "/%0.3f.png" % f_number)
