import sys
import pyHiChi as hichi
import tight_focusing_fields as sphericalPulse
import tight_focusing_show as visual
import tight_focusing_write_file as fileWriter
import hichi_primitives
import math as ma


# ------------------- initializing -------------------------------------


dirResult = "./results/"  # output directory

factor = 1
NxFull = factor*320                                # size of grid in full area
Ny = factor*256
Nz = factor*256
NxBand = 32*factor                                 # size of grid in the band
gridSize = hichi.vector3d(NxBand, Ny, Nz)          # real size of grid

# creating of spherical pulse
sphericalPulse.createSphericalPulseC(F_number_ = 1.07)

timeStep = 1*sphericalPulse.wavelength/hichi.c     # time step
maxIter = 20                                       # number of iterations to compute       

minCoords = hichi.vector3d(-20e-4, -20e-4, -20e-4) # bounds of full area
maxCoords = hichi.vector3d(20e-4, 20e-4, 20e-4)

D = 2*sphericalPulse.pulseLength                   # width of band
                                                   # to compute the task in the full area
                                                   # just set D = maxCoords.x - minCoords.x
                                                   
                                                   
# creating of usual mapping 
mapping = hichi.TightFocusingMapping(sphericalPulse.R0, sphericalPulse.pulseLength, D)

# creating of mapping with cutting at angle
#cutAngle = sphericalPulse.openingAngle + sphericalPulse.edgeSmoothingAngle
#mapping = hichi.TightFocusingMapping(sphericalPulse.R0, sphericalPulse.pulseLength, D, cutAngle)

xMin = mapping.getxMin()  # bounds if the band
xMax = mapping.getxMax()

# computing of step of grid
gridMinCoords = hichi.vector3d(xMin, minCoords.y, minCoords.z)
gridMaxCoords = hichi.vector3d(xMax, maxCoords.y, maxCoords.z)
gridStep = hichi.vector3d((gridMaxCoords.x - gridMinCoords.x) / gridSize.x, \
                          (gridMaxCoords.y - gridMinCoords.y) / gridSize.y, \
                          (gridMaxCoords.z - gridMinCoords.z) / gridSize.z)

# greating of grid for PSATDGridMapping
# 'ifCut' is a parameter to cut secondary pulses or not
ifCut = False
grid = hichi.PSATDGridMapping(gridSize, timeStep, gridMinCoords, gridStep, mapping, ifCut)
   
# creating of field solver PSATD for existing grid
fieldSolver = hichi.PSATDWithPoisson(grid)

    
def initialize():
    # mapping for tight focusing has a external parameter - time
    # so we need to advance time every iteration and to set it null in the beginning
    mapping.setTime(0.0)  
    
    # setting of start conditions
    sphericalPulse.setField(grid)


# ------------------- running -----------------------------------------------


# function to update field
def updateFields():
    global mapping
    mapping.advanceTime(timeStep)  # mapping for tight focusing has a external parameter - time
                                   # so we need to advance time every iteration and to set it null in the beginning
    fieldSolver.updateFields()     # doing one iteration of PSATD


# ----------- run and show animation (size of grid should be not large) ------
#initialize()
#visual.initVisual(minCoords, maxCoords, NxFull*2, Ny*2)
#visual.animate(grid, updateFields, maxIter=maxIter)

# ----------- run and save pictures for every iteration ----------------------
initialize()
visual.initVisual(minCoords, maxCoords, NxFull*2, Ny*2)
visual.savePictures(grid, updateFields, maxIter=maxIter, dirResult = dirResult)

# ----------- run and save results 2d in .csv files --------------------------
#initialize()
#hichi_primitives.createDir("./2d_main/")
#fileWriter.writeXOY(grid, updateFields, minCoords, maxCoords, NxFull*2, Ny*2, maxIter = maxIter, dumpIter = 1, dirResult = "./2d_main/")

# ----------- run and save results 1d in .csv files --------------------------
#initialize()
#hichi_primitives.createDir("./1d_main/")
#fileWriter.writeOX(grid, updateFields, minCoords, maxCoords, NxFull*2, maxIter = maxIter, dumpIter = 1, dirResult = "./1d_main/")




