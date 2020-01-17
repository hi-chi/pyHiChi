import sys
sys.path.append("../../build/src/pyHiChi/Release")
import pyHiChi as hichi
import tight_focusing_fields as sphericalPulse
import tight_focusing_show as visual
import tight_focusing_write_file as fileWriter
import hichi_primitives
import math as ma

factor = 1
NxFull = factor*320
Ny = factor*128
Nz = factor*128
NxBand = 128*factor #factor*32
gridSize = hichi.vector3d(NxBand, Ny, Nz)

dirResult = "./results/"

timeStep = 0.5*sphericalPulse.wavelength/hichi.c
maxIter = 64

minCoords = hichi.vector3d(-20e-4, -20e-4, -20e-4)
maxCoords = hichi.vector3d(20e-4, 20e-4, 20e-4)

#D = 2*sphericalPulse.pulseLength
D = maxCoords.x - minCoords.x

#cutAngle = sphericalPulse.openingAngle + sphericalPulse.edgeSmoothingAngle
#mapping = hichi.TightFocusingMapping(sphericalPulse.R0, sphericalPulse.pulseLength, D, cutAngle)
mapping = hichi.TightFocusingMapping(sphericalPulse.R0, sphericalPulse.pulseLength, D)

xMin = mapping.getxMin()
xMax = mapping.getxMax()

gridMinCoords = hichi.vector3d(xMin, minCoords.y, minCoords.z)
gridMaxCoords = hichi.vector3d(xMax, maxCoords.y, maxCoords.z)

gridStep = hichi.vector3d((gridMaxCoords.x - gridMinCoords.x) / gridSize.x, \
                          (gridMaxCoords.y - gridMinCoords.y) / gridSize.y, \
                          (gridMaxCoords.z - gridMinCoords.z) / gridSize.z)

ifCut = False
grid = hichi.PSATDGridMapping(gridSize, timeStep, gridMinCoords, gridStep, mapping, ifCut)

sphericalPulse.setField(grid)

fieldSolver = hichi.PSATDWithPoisson(grid)

def updateFields():
    global mapping
    mapping.advanceTime(timeStep)
    fieldSolver.updateFields()
    
visual.initVisual(minCoords, maxCoords, NxFull*5, Ny*5)
visual.animate(grid, updateFields, maxIter=maxIter)
#visual.savePictures(grid, updateFields, maxIter=maxIter, dirResult = "./pictures/")


#hichi_primitives.createDir(dirResult)

#fileWriter.write(grid, updateFields, minCoords, maxCoords, 300, 300, maxIter=160, dumpIter = 20)

#fileWriter.writeOX(grid, updateFields, minCoords, maxCoords, 64, maxIter=80, dumpIter=80, dirResult = dirResult)




