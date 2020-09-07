import sys
sys.path.append("../../python_modules")
sys.path.append("../../bin")
import pyHiChi as hichi
from tight_focusing_fields import SphericalPulseC
import numpy as np
import plotly.graph_objects as go


sphericalPulse = SphericalPulseC(f_number = 0.3,
                                 R0 = 16,
                                 pulselength = 2.0,
                                 phase = 0,
                                 edgeSmoothingAngle = 0.0
                                )

D = 2*sphericalPulse.pulselength

minFullCoords = hichi.vector3d(-20, -20, -20)
maxFullCoords = hichi.vector3d(20, 20, 20)

NxFull = 320
NxBand = D / (maxFullCoords.x - minFullCoords.x) * NxFull
Ny = 320
Nz = 320

mapping = hichi.TightFocusingMapping(sphericalPulse.R0, sphericalPulse.pulselength, D)

xMin = mapping.getxMin()
xMax = mapping.getxMax()
minBandCoords = hichi.vector3d(xMin, minFullCoords.y, minFullCoords.z)
maxBandCoords = hichi.vector3d(xMax, maxFullCoords.y, maxFullCoords.z)

gridSize = hichi.vector3d(NxBand, Ny, Nz)
gridStep = (maxBandCoords - minBandCoords) / gridSize
timeStep = 1.0/hichi.c

grid = hichi.PSATDGrid(gridSize, timeStep, minBandCoords, gridStep)
gridMapping = hichi.PSATDGridMapping(grid)
gridMapping.setMapping(mapping)

sphericalPulse.setField(gridMapping)

def get_field(grid, X, Y, Z):
    field = np.zeros(shape=X.shape)
    for i in range(X.shape[0]):
        coord = hichi.vector3d(X[i], Y[i], Z[i])
        field[i] = grid.getE(coord).norm()
    return field
    
factor = 1.0/80.0


# -------- full graph ----------

X, Y, Z = np.mgrid[minFullCoords.x:maxFullCoords.x:10j,#(int(NxFull*factor)*1j), \
                   minFullCoords.y:maxFullCoords.y:10j,#(int(Ny*factor)*1j),
                   minFullCoords.z:maxFullCoords.z:10j,#(int(Nz*factor)*1j)
                  ]
X = X.flatten()
Y = Y.flatten()
Z = Z.flatten()

values = get_field(gridMapping, X, Y, Z)
#print(values.shape)

fig = go.Figure(data=go.Isosurface(
    x=X,
    y=Y,
    z=Z,
    value=values,
    isomin=0.001,
    isomax=0.5,
    #surface_count=5, # number of isosurfaces, 2 by default: only min and max
    #colorbar_nticks=5, # colorbar ticks correspond to isosurface values
    #caps=dict(x_show=False, y_show=False)
    ))
fig.show()


# import plotly.graph_objects as go
# import numpy as np

# X, Y, Z = np.mgrid[-5:5:8j, -5:5:8j, -5:5:8j]

# # ellipsoid
# values = X * X * 0.5 + Y * Y + Z * Z * 2

# fig = go.Figure(data=go.Isosurface(
    # x=X.flatten(),
    # y=Y.flatten(),
    # z=Z.flatten(),
    # value=values.flatten(),
    # isomin=10,
    # isomax=40,
    # #caps=dict(x_show=False, y_show=False)
    # ))
# fig.show()
