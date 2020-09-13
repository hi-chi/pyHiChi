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
spherical_pulse = SphericalPulseC(f_number=0.3,
                                 R0=16,
                                 pulselength=2.0,
                                 phase=0,
                                 edge_smoothing_angle=0.1
                                )

# the computational area (coordinates of the opposite corners of the parallelepiped)
min_coords = hichi.Vector3d(-20, -20, -20)
max_coords = hichi.Vector3d(20, 20, 20)

# the field size for the entire original computational domain
FACTOR = 1.0  # the coefficient to adjust the field size
NX_FULL = int(320*FACTOR)
NY = int(256*FACTOR)
NZ = int(256*FACTOR)
full_grid_size = hichi.Vector3d(NX_FULL, NY, NZ)

# the spatial step
grid_step = (max_coords - min_coords) / full_grid_size

# the time step which is equal to R0 in the universal system of units
time_step = spherical_pulse.R0/hichi.c

# defining the minimal and the maximal width of the band (Dmin, Dmax)
Rmax = spherical_pulse.R0 + 0.5*spherical_pulse.pulselength
Rmin = spherical_pulse.R0 - 0.5*spherical_pulse.pulselength
angle = spherical_pulse.opening_angle + spherical_pulse.edge_smoothing_angle
Dmin = -Rmin*np.cos(angle) + (Rmax**2 - (Rmin*np.sin(angle))**2)**0.5
Dmax = (max_coords.x - min_coords.x)*0.5

# defining of the array with different band widths D
# the x axis of the final graph
nx_arr = np.arange(int(Dmin/grid_step.x) + 1, int(Dmax/grid_step.x) + 1, 8*FACTOR)
D_arr = nx_arr * grid_step.x


# ---------------------- run functions ----------------------  


# the function to read data from the field
def get_fields(field):
    x = np.arange(min_coords.x, max_coords.x, (max_coords.x - min_coords.x)/NX_FULL)
    y = np.arange(min_coords.y, max_coords.y, (max_coords.y - min_coords.y)/NY)
    res = np.zeros(shape=(NX_FULL, NY))
    for ix in range(NX_FULL):
        for iy in range(NY):
            coord_xy = hichi.Vector3d(x[ix], y[iy], 0.0)
            E = field.get_E(coord_xy)
            res[ix, iy] = E.norm()
    return res


# the function to run the simulation
def run(nx):

    # the band width
    D = nx * grid_step.x

    # the mapping
    mapping = hichi.TightFocusingMapping(spherical_pulse.R0,
                                         spherical_pulse.pulselength,
                                         D
                                        )
    
    # the bounds of the band
    x_min = mapping.get_min_coord()
    x_max = mapping.get_max_coord()   
    band_min_coords = hichi.Vector3d(x_min, min_coords.y, min_coords.z)
    band_max_coords = hichi.Vector3d(x_max, max_coords.y, max_coords.z)
    
    # the field size for the periodic space (that is the real field size)
    grid_size = hichi.Vector3d(nx, NY, NZ)
    
    # the creation of the computational field with the approciate mapping
    field = hichi.PSATDField(grid_size, band_min_coords, grid_step, time_step)
    field = field.apply_mapping(mapping)
    
    # the field initialisation
    spherical_pulse.set_field(field)
    
    # the correction of the start conditions using the Poisson's equation
    field.convert_fields_poisson_equation()
    
    # the performing one iteration of the field solver
    field.update_fields()

    return get_fields(field)


# the function to compute error of the scheme
# returns the maximal difference between full and reduced calculations
def compute_error(res_band, res_full):
    m = 0.0
    for iy in range(NY):
        for ix in range(NX_FULL):
            diff = abs(res_band[ix, iy] - res_full[ix, iy])
            m = diff if m < diff else m
    return m


# ---------------------- run -----------------------

res_error = []  # the y axis of the final graph

res_full = run(NX_FULL)  # the full simulation result

for nx in nx_arr:
    print("%d/%d\r" % (np.where(nx_arr == nx)[0], len(nx_arr)), end="")
    res_band = run(nx)  # the band simulation result
    res_error.append(compute_error(res_band, res_full))
    

# ---------------------- plot graph------------------

import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams.update({"font.size" : 17})

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(D_arr/spherical_pulse.pulselength, res_error/res_full.max()*100, "->b")
ax.set_xlabel('$D/L$')
ax.set_ylabel('error, %')
ax.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator([Dmin/spherical_pulse.pulselength]))
ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(2))
ax.yaxis.set_minor_locator(matplotlib.ticker.FixedLocator([0]))
ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))
ax.grid(which='major', linestyle='--')
ax.grid(which='minor', linestyle='-')
fig.tight_layout()
plt.savefig(DIR_RESULT + "/error_D.png")

