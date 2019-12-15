import sys
sys.path.append("../../build/src/pyHiChi/Release")
import pyHiChi as hichi
import tight_focusing_fields as sphericalPulse
import tight_focusing_show as visual

factor = 1
gridSize = hichi.vector3d(factor*64, factor*64, factor*32)

widthCoeff = 1.2
coordMax = widthCoeff*sphericalPulse.R0
minCoords = hichi.vector3d(-coordMax, -coordMax, -coordMax)
maxCoords = hichi.vector3d(coordMax, coordMax, coordMax)

D = 2*sphericalPulse.pulseLength

xMin = sphericalPulse.R0 - 0.5*D
xMax = sphericalPulse.R0 + 0.5*D

gridStep = hichi.vector3d((maxCoords.x - minCoords.x) / gridSize.x, \
                          (maxCoords.y - minCoords.y) / gridSize.y, \
                          (maxCoords.z - minCoords.z) / gridSize.z)

timeStep = 0.2*sphericalPulse.wavelength/hichi.c

#mapping = hichi.Mapping()

grid = hichi.PSATDGrid(gridSize, timeStep, minCoords, gridStep)
#grid.setMapping(mapping)

fieldFuncs = sphericalPulse.getFieldFuncs()
grid.setE(fieldFuncs[0][0], fieldFuncs[0][1], fieldFuncs[0][2])
grid.setB(fieldFuncs[1][0], fieldFuncs[1][1], fieldFuncs[1][2])

fieldSolver = hichi.PSATD(grid)

def updateFields():
    fieldSolver.updateFields()
    
visual.initVisual(minCoords, maxCoords)
visual.show(grid, updateFields)





