import sys
sys.path.append("../")
import pyHiChi as hichi
import tight_focusing_fields as sphericalPulse
import tight_focusing_show as visual
import tight_focusing_write_file as fileWriter
import hichi_primitives
import math as ma
import time

# ------------------- initializing -------------------------------------


factor = 0.5
NxFull = int(factor*320)                           # size of grid in full area
Ny = int(factor*256)
Nz = int(factor*256)

# creating of spherical pulse
sphericalPulse.createSphericalPulseC(F_number_ = 0.3, R0_ = 16, pulselength_ = 2.0, phase_ = 0, edgeSmoothingAngle_ = 0.1)

timeStep = sphericalPulse.getDtCGS(1.0)            # time step in CGS
maxIter = 32                                       # number of iterations to compute       

minCoords = hichi.vector3d(-20, -20, -20)          # bounds of full area
maxCoords = hichi.vector3d(20, 20, 20)


D_div_L_arr = range(2, 22, 2)


for D_div_L in D_div_L_arr:

    D = D_div_L*sphericalPulse.pulselength             # width of band
    print("D = %d L" % D_div_L)
    sys.stdout.flush()    
    
    NxBand = int(D*NxFull/(maxCoords.x-minCoords.x))   # size of grid in the band
    gridSize = hichi.vector3d(NxBand, Ny, Nz)          # real size of grid
    
    START_INIT = time.time()
    
    mapping = hichi.TightFocusingMapping(sphericalPulse.R0, sphericalPulse.pulselength, D)
    
    xMin = mapping.getxMin()
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
    
    fieldSolver = hichi.PSATD(grid)
    
    def initialize():
        mapping.setTime(0.0)  
        sphericalPulse.setField(grid)
        fieldSolver.convertFieldsPoissonEquation()
    
    def update():
        mapping.advanceTime(timeStep)   
        fieldSolver.updateFields()

    initialize()
    
    END_INIT = time.time()
    print("Python: time of init is ", (END_INIT - START_INIT)*1000, " ms")
    sys.stdout.flush()
    
    for iter in range(maxIter):
        START_ITER = time.time()
        update()
        END_ITER = time.time()
        print("Python: time of %d iter is " % iter, (END_ITER - START_ITER)*1000, " ms")
        sys.stdout.flush()

