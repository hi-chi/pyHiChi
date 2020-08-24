import sys
sys.path.append("../../python_modules")
sys.path.append("../../bin")

import pyHiChi as hichi
import numpy as np
from tight_focusing_fields import SphericalPulseC

# the output directory
DIR_RESULT = "./"

# the creation of the spherical pulse
# f_number=0.3 (opening angle = 1 rad)
sphericalPulse = SphericalPulseC(f_number=0.3,
                                 R0=16,
                                 pulselength=2.0,
                                 phase=0,
                                 edgeSmoothingAngle=0.1
                                )

# the computational area (coordinates of the opposite corners of the parallelepiped)
minCoords = hichi.vector3d(-20, -20, -20)
maxCoords = hichi.vector3d(20, 20, 20)

# the grid size for the entire original computational domain
factor = 1.0  # the coefficient to adjust the grid size
NxFull = int(320*factor)
Ny = int(256*factor)
Nz = int(256*factor)
fullGridSize = hichi.vector3d(NxFull, Ny, Nz)

# the spatial step
gridStep = (maxCoords - minCoords) / fullGridSize

# the time step which is equal to R0 in the universal system of units
timeStep = sphericalPulse.R0/hichi.c

# defining the minimal and the maximal width of the band (Dmin, Dmax)
Rmax = sphericalPulse.R0 + 0.5*sphericalPulse.pulselength
Rmin = sphericalPulse.R0 - 0.5*sphericalPulse.pulselength
angle = sphericalPulse.openingAngle + sphericalPulse.edgeSmoothingAngle
Dmin = -Rmin*np.cos(angle) + (Rmax**2 - (Rmin*np.sin(angle))**2)**0.5
Dmax = (maxCoords.x - minCoords.x)*0.5

# defining of the array with different band widths D
# the x axis of the final graph
NxArr = np.arange(int(Dmin/gridStep.x) + 1, int(Dmax/gridStep.x) + 1, 8*factor)
DArr = NxArr * gridStep.x


# ---------------------- run functions ----------------------  


# the function to read data from the grid
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


# the function to run the simulation
def run(Nx):

    # the band width
    D = Nx * gridStep.x

    # the mapping
    mapping = hichi.TightFocusingMapping(sphericalPulse.R0,
                                         sphericalPulse.pulselength,
                                         D)
    
    # the bounds of the band
    xMin = mapping.getxMin()
    xMax = mapping.getxMax()   
    gridMinCoords = hichi.vector3d(xMin, minCoords.y, minCoords.z)
    gridMaxCoords = hichi.vector3d(xMax, maxCoords.y, maxCoords.z)
    
    # the grid size for the periodic space (that is the real grid size)
    gridSize = hichi.vector3d(Nx, Ny, Nz)
    
    # the creation of the computational grid with the approciate mapping
    grid = hichi.PSATDGridMapping(gridSize, timeStep, gridMinCoords, gridStep)
    grid.setMapping(mapping)
    
    # the field solver (PSATD)
    fieldSolver = hichi.PSATD(grid)
    
    # the grid initialisation
    sphericalPulse.setField(grid)
    
    # the correction of the start conditions using the Poisson's equation
    fieldSolver.convertFieldsPoissonEquation()
    
    # the performing one iteration of the field solver
    fieldSolver.updateFields()

    return getFields(grid)
    

# the function to compute error of the scheme
# returns the maximal difference between full and reduced calculations
def computeError(resBand, resFull):
    m = 0.0
    for iy in range(Ny):
        for ix in range(NxFull):
            diff = abs(resBand[ix, iy] - resFull[ix, iy])
            m = diff if m < diff else m
    return m


# ---------------------- run -----------------------

resError = []  # the y axis of the final graph

resFull = run(NxFull)  # the full simulation result

for Nx in NxArr:
    print("%d/%d\r" % (np.where(NxArr == Nx)[0], len(NxArr)), end="")
    resBand = run(Nx)  # the band simulation result
    resError.append(computeError(resBand, resFull))
    

# ---------------------- plot graph------------------

import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams.update({"font.size" : 17})

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(DArr/sphericalPulse.pulselength, resError/resFull.max()*100, "->b")
ax.set_xlabel('$D/L$')
ax.set_ylabel('error, %')
ax.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator([Dmin/sphericalPulse.pulselength]))
ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(2))
ax.yaxis.set_minor_locator(matplotlib.ticker.FixedLocator([0]))
ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))
ax.grid(which='major', linestyle='--')
ax.grid(which='minor', linestyle='-')
fig.tight_layout()
plt.savefig(DIR_RESULT + "/error_D.png")

