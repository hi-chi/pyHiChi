import sys
sys.path.append("../../python_modules")
sys.path.append("../../bin")
import pyHiChi as hichi
from tight_focusing_fields import SphericalPulseC, SphericalPulsePython


# ------------------- initializing -------------------------------------


FACTOR = 0.5
NX_FULL = int(FACTOR*320)                           # size of field in full area
NY = int(FACTOR*256)
NZ = int(FACTOR*256)
NX_BAND = int(56*FACTOR)                            # size of field in the band
grid_size = hichi.Vector3d(NX_BAND, NY, NZ)          # real size of field

# creating of spherical pulse
spherical_pulse = SphericalPulsePython(f_number = 0.3,
                                       R0 = 16,
                                       pulselength = 2.0,
                                       phase = 0,
                                       edge_smoothing_angle = 0.1
                                      )

time_step = 1.0/hichi.c                             # time step in CGS system of units
max_iter = 32                                       # number of iterations to compute       

min_coords = hichi.Vector3d(-20, -20, -20)          # bounds of full area
max_coords = hichi.Vector3d(20, 20, 20)

D = 3.5*spherical_pulse.pulselength                 # band width

# to compute the task in the full area just set
# D = max_coords.x - min_coords.x
                                                   
                                                   
# creating of mapping
mapping = hichi.TightFocusingMapping(spherical_pulse.R0, spherical_pulse.pulselength, D)

# not to cut secondary pulses
# mapping.if_perform_inverse_mapping(False)

x_min = mapping.get_min_coord()  # bounds of the band
x_max = mapping.get_max_coord()

# computing of step of field
band_min_coords = hichi.Vector3d(x_min, min_coords.y, min_coords.z)
band_max_coords = hichi.Vector3d(x_max, max_coords.y, max_coords.z)
grid_step = (band_max_coords - band_min_coords) / grid_size

# creating of field for PSATDGrid
field = hichi.PSATDPoissonField(grid_size, band_min_coords, grid_step, time_step)
field = field.apply_mapping(mapping)


def initialize():
    # setting of start conditions
    spherical_pulse.set_field(field)


def update_fields():
    # doing one iteration of PSATD
    field.update_fields()


# ----------- run and show animation (size of field should be not large) ------

from hichi_visualisation import *
from hichi_primitives import Axis, Plane, Field

visual = Visual(field, min_coords, max_coords, dpi=500, fontsize=17)

def animate_plane(visual, n_iter):
    visual.animate_plane(shape=(NX_FULL*2, NY*2), func_update=update_fields, n_iter=n_iter,
                          plane=Plane.XOY, last_coordinate_value=0.0,
                          field=Field.E, norm=True,
                          value_limits=(0.0, 0.5),
                          interval=1
                         )
                            
def animate_axis(visual, n_iter):
    visual.animate_axis(n_points=NX_FULL*2, func_update=update_fields, n_iter=n_iter,
                         axis=Axis.X, last_coordinate_value=(0.0, 0.0),
                         field=Field.E, norm=True,
                         y_limits=(-1.0, 8.0),
                         interval=50
                        )                            
                            
initialize()
animate_plane(visual, max_iter)  # it should be the last function in the current script
#animate_axis(visual, max_iter)  # it should be the last function in the current script


# ----------- run and save pictures for every iteration ----------------------
'''
from hichi_visualisation import *
from hichi_primitives import Axis, Plane, Field, create_dir

create_dir("./pictures")
visual = Visual(field, min_coords, max_coords, "./pictures", dpi=500, fontsize=17)

def save_plane_to_image(visual, iter):
    visual.save_plane_to_image(shape=(NX_FULL*4, NY*4), plane=Plane.XOY, last_coordinate_value=0.0,
                               field=Field.E, norm=True,
                               value_limits=(0.0, 0.5),
                               name_picture="field%04d.png" % iter
                              )
                            
def save_axis_to_image(visual, iter):
    visual.save_axis_to_image(n_points=NX_FULL*4, axis=Axis.X, last_coordinate_value=(0.0, 0.0),
                              field=Field.E, norm=True,
                              y_limits=(-1.0, 8.0),
                              name_picture="field%04d.png" % iter
                             )                            
                            
initialize()
for i in range(max_iter):
    save_plane_to_image(visual, i)
    #save_axis_to_image(visual, i)
    update_fields()
'''

# ----------- run and save results in .csv files --------------------------
'''
from hichi_writing import Writer, Reader
from hichi_primitives import Axis, Plane, Field, create_dir

create_dir("./csv")
writer = Writer(field, min_coords, max_coords, "./csv")

def save_plane_to_file(writer, iter):
    writer.save_plane_to_file(shape=(NX_FULL*4, NY*4), plane=Plane.XOY, last_coordinate_value=0.0,
                              field=Field.E, norm=True,
                              name_file="field%04d.csv" % iter
                             )
                            
def save_axis_to_file(writer, iter):
    writer.save_axis_to_file(n_points=NX_FULL*4, axis=Axis.X, last_coordinate_value=(0.0, 0.0),
                             field=Field.E, norm=True,
                             name_file="field%04d.csv" % iter
                            )                            
                            
initialize()
for i in range(max_iter):
    save_plane_to_file(writer, i)
    #save_axis_to_file(writer, i)
    update_fields()


# it is possible to read crated file to numpy.array   
import numpy as np
reader = Reader("./csv")
field = reader.read_file_2d("field0000.csv")
#field = reader.read_file_1d("field0000.csv")
print(field)
'''
