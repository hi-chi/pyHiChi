import sys
sys.path.append("../../build/src/pyHiChi/Release")
import pyHiChi as hichi
import tight_focusing_fields as sphericalPulse
import tight_focusing_show as visual
import math as ma

factor = 1
gridSize = hichi.vector3d(factor*32, factor*128, factor*128)

timeStep = 0.2*sphericalPulse.wavelength/hichi.c

D = 2*sphericalPulse.pulseLength

#cutAngle = sphericalPulse.openingAngle + sphericalPulse.edgeSmoothingAngle
#mapping = hichi.TightFocusingMapping(sphericalPulse.R0, sphericalPulse.pulseLength, D, cutAngle)
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

ifCut = False
grid = hichi.PSATDGridMapping(gridSize, timeStep, gridMinCoords, gridStep, mapping, ifCut)

fieldFuncs = sphericalPulse.getFieldFuncs()
grid.setE(fieldFuncs[0][0], fieldFuncs[0][1], fieldFuncs[0][2])
grid.setB(fieldFuncs[1][0], fieldFuncs[1][1], fieldFuncs[1][2])

fieldSolver = hichi.PSATD(grid)

def updateFields():
    global mapping
    mapping.advanceTime(timeStep)
    fieldSolver.updateFields()

visual.initVisual(minCoords, maxCoords)
visual.show(grid, updateFields)





