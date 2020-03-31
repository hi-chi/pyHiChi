import sys
import pyHiChi as hichi
import tight_focusing_fields as sphericalPulse
import tight_focusing_show as visual
import tight_focusing_write_file as fileWriter
import hichi_primitives
import math as ma


# ------------------- initializing -------------------------------------


factor = 0.5
NxFull = int(factor*320)                           # size of grid in full area
Ny = int(factor*256)
Nz = int(factor*256)
NxBand = NxFull#int(56*factor)                            # size of grid in the band
gridSize = hichi.vector3d(NxBand, Ny, Nz)          # real size of grid

# creating of spherical pulse
sphericalPulse.createSphericalPulsePython(F_number_ = 0.3, R0_ = 16, pulselength_ = 2.0, phase_ = 0, edgeSmoothingAngle_ = 0.1)

timeStep = sphericalPulse.getDtCGS(1.0)            # time step in CGS
maxIter = 32                                       # number of iterations to compute       

minCoords = hichi.vector3d(-20, -20, -20)          # bounds of full area
maxCoords = hichi.vector3d(20, 20, 20)

D = 3.5*sphericalPulse.pulselength                 # width of band

# to compute the task in the full area just set
D = maxCoords.x - minCoords.x
                                                   
                                                   
# creating of mapping
mapping = hichi.TightFocusingMapping(sphericalPulse.R0, sphericalPulse.pulselength, D)

# creating of mapping with cutting at angle
#cutAngle = sphericalPulse.openingAngle + sphericalPulse.edgeSmoothingAngle
#mapping = hichi.TightFocusingMapping(sphericalPulse.R0, sphericalPulse.pulseLength, D, cutAngle)

# not to cut secondary pulses
mapping.setIfCut(False)

xMin = mapping.getxMin()  # bounds of the band
xMax = mapping.getxMax()

# computing of step of grid
gridMinCoords = hichi.vector3d(xMin, minCoords.y, minCoords.z)
gridMaxCoords = hichi.vector3d(xMax, maxCoords.y, maxCoords.z)
gridStep = hichi.vector3d((gridMaxCoords.x - gridMinCoords.x) / gridSize.x, \
                          (gridMaxCoords.y - gridMinCoords.y) / gridSize.y, \
                          (gridMaxCoords.z - gridMinCoords.z) / gridSize.z)

# creating of grid for PSATDGridMapping
grid = hichi.PSATDGridMapping(gridSize, timeStep, gridMinCoords, gridStep)
grid.setMapping(mapping)

# creating of field solver PSATD for existing grid
fieldSolver = hichi.PSATD(grid)


def initialize():
    # mapping for tight focusing has a external parameter - time
    # so we need to advance time every iteration and to set it null in the beginning
    mapping.setTime(0.0)  
    
    # setting of start conditions
    sphericalPulse.setField(grid)
    
    # consider the Poisson equation
    fieldSolver.convertFieldsPoissonEquation()


# ------------------- running -----------------------------------------------


# function to update field
def updateFields():
    mapping.advanceTime(timeStep)  # mapping for tight focusing has an external parameter time
                                   # so we need to advance time every iteration and to set it null in the beginning
    fieldSolver.updateFields()     # doing one iteration of PSATD


# ----------- run and show animation (size of grid should be not large) ------
#initialize()
#visual.initVisual(minCoords, maxCoords, NxFull*2, Ny*2)
#visual.animate(grid, updateFields, maxIter=maxIter)

# ----------- run and save pictures for every iteration ----------------------
initialize()
visual.initVisual(minCoords, maxCoords, NxFull*4, Ny*4)
visual.savePictures(grid, updateFields, maxIter=maxIter, dirResult = "./pictures/")  # DELETE OLD FILES IN "./pictures/"!!!

# ----------- run and save results 2d in .csv files --------------------------
#initialize()
#hichi_primitives.createDir("./2d_main/")  # DELETE OLD FILES IN ./2d_main/!!!
#fileWriter.writeXOY(grid, updateFields, minCoords, maxCoords, NxFull*2, Ny*2, maxIter = maxIter, dumpIter = 1, dirResult = "./2d_main/")

# ----------- run and save results 1d in .csv files --------------------------
#initialize()
#hichi_primitives.createDir("./1d_main/")  # DELETE OLD FILES IN ./1d_main/!!!
#fileWriter.writeOX(grid, updateFields, minCoords, maxCoords, NxFull*2, maxIter = maxIter, dumpIter = 1, dirResult = "./1d_main/")




