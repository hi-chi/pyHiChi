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
# 1st arg is axis around which rotation is
# 2nd arg is angle
mapping1 = hichi.RotationMapping(hichi.Axis.X, np.pi/2)
field1 = field + field.apply_mapping(mapping1)
# field.apply_mapping(mapping1) has another polarization vector
# fields have not to be destroyed

# create the second transformed field
mapping2 = hichi.RotationMapping(hichi.Axis.X, np.pi)
field2 = field + field.apply_mapping(mapping2)
# fields have to be destroyed

# create the first transformed field
# 3rd arg is propagation direction (x)
# 4th arg is an additional angle by which the polarization vector rotates around the propagation direction
mapping3 = hichi.RotationMapping(hichi.Axis.X, np.pi, hichi.Axis.X, np.pi/2)
field3 = field + field.apply_mapping(mapping3)
# fields have not to be destroyed

# create the first transformed field
mapping4 = hichi.RotationMapping(hichi.Axis.X, np.pi/2, hichi.Axis.X, np.pi/2)
field4 = field + field.apply_mapping(mapping4)
# fields have to be destroyed


def update_data():
    field.update_fields()
    
    
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

# show field1
im1 = axes[0, 0].imshow(get_fields(field1), cmap='RdBu', interpolation='none', vmax=1.1,
    extent=(min_coords.x, max_coords.x, min_coords.y, max_coords.y), animated = True)
fig.colorbar(im1, ax=axes[0, 0])
axes[0, 0].set_xlabel("x")
axes[0, 0].set_ylabel("y")
axes[0, 0].set_title("field1")

# show field2
im2 = axes[0, 1].imshow(get_fields(field2), cmap='RdBu', interpolation='none', vmax=1.1,
    extent=(min_coords.x, max_coords.x, min_coords.y, max_coords.y), animated = True)
fig.colorbar(im2, ax=axes[0, 1])
axes[0, 1].set_xlabel("x")
axes[0, 1].set_ylabel("y")
axes[0, 1].set_title("field2")

# show field3
im3 = axes[1, 0].imshow(get_fields(field3), cmap='RdBu', interpolation='none', vmax=1.1,
    extent=(min_coords.x, max_coords.x, min_coords.y, max_coords.y), animated = True)
fig.colorbar(im3, ax=axes[1, 0])
axes[1, 0].set_xlabel("x")
axes[1, 0].set_ylabel("y")
axes[1, 0].set_title("field3")

# show field4
im4 = axes[1, 1].imshow(get_fields(field4), cmap='RdBu', interpolation='none', vmax=1.1,
    extent=(min_coords.x, max_coords.x, min_coords.y, max_coords.y), animated = True)
fig.colorbar(im2, ax=axes[1, 1])
axes[1, 1].set_xlabel("x")
axes[1, 1].set_ylabel("y")
axes[1, 1].set_title("field4")


def update_fig(*args):
    update_data()
    im1.set_array(get_fields(field1))
    im2.set_array(get_fields(field2))
    im3.set_array(get_fields(field3))
    im4.set_array(get_fields(field4))
    return im1, im2, im3, im4,
    
ani = animation.FuncAnimation(fig, update_fig, interval=10, blit=True)

fig.tight_layout()

plt.show()
