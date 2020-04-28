import sys
sys.path.append("../")
sys.path.append("../../python_modules")
import pyHiChi as hichi
import tight_focusing_fields as sphericalPulse
import math as ma
import time

# ------------------- initializing -------------------------------------

START_SCRIPT = time.time()

factor = 3.5
NxFull = int(factor*320)                           # size of grid in full area
Ny = int(factor*256)
Nz = int(factor*256)

# creating of spherical pulse
sphericalPulse.createSphericalPulseC(F_number_ = 0.3, R0_ = 16, pulselength_ = 2.0, phase_ = 0, edgeSmoothingAngle_ = 0.1)

timeStep = sphericalPulse.getDtCGS(8.0)            # time step in CGS
maxIter = 4                                        # number of iterations to compute       

minCoords = hichi.vector3d(-20, -20, -20)          # bounds of full area
maxCoords = hichi.vector3d(20, 20, 20)


D_div_L_arr = range(2, 22, 2)


for D_div_L in D_div_L_arr:

    START_D = time.time()

    D = D_div_L*sphericalPulse.pulselength             # width of band
    print("D = %d L" % D_div_L)
    sys.stdout.flush()    
    
    NxBand = int(D*NxFull/(maxCoords.x-minCoords.x))   # size of grid in the band
    gridSize = hichi.vector3d(NxBand, Ny, Nz)          # real size of grid
    
    mapping = hichi.TightFocusingMapping(sphericalPulse.R0, sphericalPulse.pulselength, D)
    
    xMin = mapping.getxMin()
    xMax = mapping.getxMax()
    
    # computing of step of grid
    gridMinCoords = hichi.vector3d(xMin, minCoords.y, minCoords.z)
    gridMaxCoords = hichi.vector3d(xMax, maxCoords.y, maxCoords.z)
    gridStep = (gridMaxCoords - gridMinCoords) / gridSize
    
    # creating of grid for PSATDGridMapping
    grid = hichi.PSATDGridMapping(gridSize, timeStep, gridMinCoords, gridStep)
    grid.setMapping(mapping)
    
    fieldSolver = hichi.PSATD(grid)
    
    def initialize():
        START_INIT = time.time()
        sphericalPulse.setField(grid)
        END_INIT = time.time()
        print("Python: time of init is %f ms" % ((END_INIT - START_INIT)*1000))
    
    def update(): 
        fieldSolver.updateFields()

    initialize()
    sys.stdout.flush()
    
    for iter in range(maxIter):
        START_ITER = time.time()
        update()
        END_ITER = time.time()
        print("Python: time of %d iter is %f ms" % (iter, (END_ITER - START_ITER)*1000))
        sys.stdout.flush()
    
    END_D = time.time()
    print("Python: time of al is %f ms" % (END_D - START_D)*1000)

END_SCRIPT = time.time()
print("Python: time of script is %f, min" % ((END_SCRIPT - START_SCRIPT)/60))
