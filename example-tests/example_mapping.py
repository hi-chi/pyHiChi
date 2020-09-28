import sys
sys.path.append("../bin")
import pyHiChi as hichi
import numpy as np


def field_value(x, y, z):
    return np.exp(-x**2-y**2-z**2)*np.sin(5*x)
    
def null_value(x, y, z):
    return 0.0
    
min_coords = hichi.Vector3d(-5, -5, -5)
max_coords = hichi.Vector3d(5, 5, 5)

grid_size = hichi.Vector3d(64, 64, 64)
grid_step = (max_coords - min_coords) / grid_size  
time_step = 0.1/hichi.c

field = hichi.PSATDField(grid_size, min_coords, grid_step, time_step)
field.set_E(null_value, field_value, null_value)
field.set_B(null_value, null_value, field_value)
field.convert_fields_poisson_equation()


# create the first transformed field
# rotation
def getRotatedField(field):
    mapping = hichi.RotationMapping(hichi.Axis.Z, -np.pi/6)  # local object
    return field.apply_mapping(mapping)  # mapping is inaccessible but not destroyed

field1 = getRotatedField(field)  


# create some other mappings  
angle = np.pi/3
rotation_mapping_z_axis = hichi.RotationMapping(hichi.Axis.Z, angle)

shift = hichi.Vector3d(-2.0, -2.0, 0.0)
shift_mapping = hichi.ShiftMapping(shift)

scaleCoef = 1.5
scaleMapping = hichi.ScaleMapping(hichi.Axis.X, scaleCoef)


# transform a point (rotation -> shift)
point = hichi.Vector3d(-1.0, -1.0, 0.0)
point2 = shift_mapping.get_inverse_coords(rotation_mapping_z_axis.get_inverse_coords(point))
print("Direct mapping: point (%0.2f, %0.2f, %0.2f) -> (%0.2f, %0.2f, %0.2f)" % \
       (point.x, point.y, point.z, point2.x, point2.y, point2.z))
point3 = shift_mapping.get_direct_coords(rotation_mapping_z_axis.get_direct_coords(point2))
print("Inverse mapping: point (%0.2f, %0.2f, %0.2f) -> (%0.2f, %0.2f, %0.2f)" % \
       (point2.x, point2.y, point2.z, point3.x, point3.y, point3.z))


# create the second transformed field
# scale -> rotation -> shift
field2_tmp1 = field.apply_mapping(scaleMapping)
field2_tmp2 = field2_tmp1.apply_mapping(rotation_mapping_z_axis)
field2 = field2_tmp2.apply_mapping(shift_mapping)


# create the third transformed field
# shift -> rotation -> scale
field3 = field.apply_mapping(shift_mapping)\
              .apply_mapping(rotation_mapping_z_axis)\
              .apply_mapping(scaleMapping)
# intermediate fields are inaccesible but not destroyed

def update_data():
    field.update_fields()
    # == field1.update_fields()
    # == field_2.update_fields()
    # == field3.update_fields()
    
    
# --------- show -------------


import matplotlib.pyplot as plt
import matplotlib.animation as animation

N = 128
x = np.linspace(min_coords.x, max_coords.x, N)
y = np.linspace(min_coords.y, max_coords.y, N)

def get_fields(field):
    global x, y, N
    res = np.zeros(shape=(N,N))
    for ix in range(N):
        for iy in range(N):
            coord_xy = hichi.Vector3d(x[ix], y[iy], 0.0)
            res[N - iy - 1, ix] = field.get_E(coord_xy).norm()
    return res

fig, axes = plt.subplots(ncols=2, nrows=2)

# show setted field
im1 = axes[0, 0].imshow(get_fields(field), cmap='RdBu', interpolation='none',
    extent=(min_coords.x, max_coords.x, min_coords.y, max_coords.y), animated = True)
fig.colorbar(im1, ax=axes[0, 0])
axes[0, 0].set_xlabel("x")
axes[0, 0].set_ylabel("y")
axes[0, 0].set_title("Original field")


# show transformed field (rotation)
im2 = axes[0, 1].imshow(get_fields(field1), cmap='RdBu', interpolation='none',
    extent=(min_coords.x, max_coords.x, min_coords.y, max_coords.y), animated = True)
fig.colorbar(im2, ax=axes[0, 1])
axes[0, 1].set_xlabel("x")
axes[0, 1].set_ylabel("y")
axes[0, 1].set_title("Transformed field \n rotation")


# show transformed field (scale -> rotation -> shift)
im3 = axes[1, 0].imshow(get_fields(field2), cmap='RdBu', interpolation='none',
    extent=(min_coords.x, max_coords.x, min_coords.y, max_coords.y), animated = True)
fig.colorbar(im3, ax=axes[1, 0])
axes[1, 0].set_xlabel("x")
axes[1, 0].set_ylabel("y")
axes[1, 0].set_title("Transformed field \n scale -> rotation -> shift")


# show transformed field (shift + rotation)
im4 = axes[1, 1].imshow(get_fields(field3), cmap='RdBu', interpolation='none',
    extent=(min_coords.x, max_coords.x, min_coords.y, max_coords.y), animated = True)
fig.colorbar(im4, ax=axes[1, 1])
axes[1, 1].set_xlabel("x")
axes[1, 1].set_ylabel("y")
axes[1, 1].set_title("Transformed field \n shift ->rotation -> scale")


def update_fig(*args):
    update_data()
    im1.set_array(get_fields(field))
    im2.set_array(get_fields(field1))
    im3.set_array(get_fields(field2))
    im4.set_array(get_fields(field3))
    return im1, im2, im3, im4,
    
ani = animation.FuncAnimation(fig, update_fig, interval=10, blit=True)

fig.tight_layout()

plt.show()
