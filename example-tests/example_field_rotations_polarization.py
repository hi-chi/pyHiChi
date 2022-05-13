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


field0 = hichi.PSATDField(grid_size, min_coords, grid_step, time_step)
# polarization vector is p0 = (0, 1, 0)
field0.set_E(null_value, field_value, null_value)
field0.set_B(null_value, null_value, field_value)
field0.convert_fields_poisson_equation()


# create the first transformed field0
# 1st arg is the axis around which the rotation is
# 2nd arg is the rotation angle
mapping1 = hichi.RotationMapping(hichi.Axis.X, np.pi/2)
field1 = field0.apply_mapping(mapping1)
# field1 has another polarization vector
# p1 = (0, 0, 1)
field1_sum = field0 + field1
# field1_sum has not to be zero

# create the second transformed field0
mapping2 = hichi.RotationMapping(hichi.Axis.X, np.pi)
field2 = field0.apply_mapping(mapping2)
# p2 = (0, -1, 0)
field2_sum = field0 + field2
# field2_sum has to be zero

# create the first transformed field0
# 3rd arg is the propagation direction (x)
# 4th arg is the additional angle by which the polarization vector rotates around the propagation direction
mapping3 = hichi.RotationMapping(hichi.Axis.X, np.pi, hichi.Axis.X, np.pi/2)
field3 = field0.apply_mapping(mapping3)
# p3 = (0, 0, -1)
# p0 = (0, 1, 0) were rotated by 3pi/2
field3_sum = field0 + field3
# field3_sum has not to be zero

# create the first transformed field0
mapping4 = hichi.RotationMapping(hichi.Axis.X, np.pi/2, hichi.Axis.X, np.pi/2)
# p4 = (0, -1, 0)
field4 = field0.apply_mapping(mapping4)
field4_sum = field0 + field4
# field4_sum has to be zero


def update_data():
    field0.update_fields()
    
    
# --------- show -------------


import matplotlib.pyplot as plt
import matplotlib.animation as animation

N = 128
x = np.linspace(min_coords.x, max_coords.x, N)
y = np.linspace(min_coords.y, max_coords.y, N)

def get_fields(field, func):
    global x, y, N
    res = np.zeros(shape=(N,N))
    for ix in range(N):
        for iy in range(N):
            coord_xy = hichi.Vector3d(x[ix], y[iy], 0.0)
            res[N - iy - 1, ix] = func(field, coord_xy)
    return res


# ------ components of each field0 --------

for ind, field in enumerate([field0, field1, field2, field3, field4]):
    def ex_func(field, coord): return field.get_Ex(coord)
    def ey_func(field, coord): return field.get_Ey(coord)
    def ez_func(field, coord): return field.get_Ez(coord)
    def bx_func(field, coord): return field.get_Bx(coord)
    def by_func(field, coord): return field.get_By(coord)
    def bz_func(field, coord): return field.get_Bz(coord)
    arr_funcs = [[ex_func, ey_func, ez_func], [bx_func, by_func, bz_func]]
    
    fig, axes = plt.subplots(nrows=2, ncols=3)
    for i, fn in enumerate(['E', 'B']):
        for j, cn in enumerate(['x', 'y', 'z']):
            im = axes[i, j].imshow(get_fields(field, arr_funcs[i][j]), cmap='RdBu', interpolation='none',
                vmin=-1.0, vmax=1.0,
                extent=(min_coords.x, max_coords.x, min_coords.y, max_coords.y), animated = True)
            fig.colorbar(im, ax=axes[i, j])
            axes[i, j].set_xlabel("x")
            axes[i, j].set_ylabel("y")
            axes[i, j].set_title(fn + cn)
    fig.suptitle("field%d" % (ind))
    fig.tight_layout()
    plt.show()
    plt.close(fig)


# ------ animation of summed fields --------

fig, axes = plt.subplots(ncols=2, nrows=2)

def norm_func(field, coord):
    return field.get_E(coord).norm()

# show field1
im1 = axes[0, 0].imshow(get_fields(field1_sum, norm_func), cmap='RdBu', interpolation='none',
    vmin=0.0, vmax=1.1,
    extent=(min_coords.x, max_coords.x, min_coords.y, max_coords.y), animated = True)
fig.colorbar(im1, ax=axes[0, 0])
axes[0, 0].set_xlabel("x")
axes[0, 0].set_ylabel("y")
axes[0, 0].set_title("field0 + field1")

# show field2
im2 = axes[0, 1].imshow(get_fields(field2_sum, norm_func), cmap='RdBu', interpolation='none',
    vmin=0.0, vmax=1.1,
    extent=(min_coords.x, max_coords.x, min_coords.y, max_coords.y), animated = True)
fig.colorbar(im2, ax=axes[0, 1])
axes[0, 1].set_xlabel("x")
axes[0, 1].set_ylabel("y")
axes[0, 1].set_title("field0 + field2")

# show field3
im3 = axes[1, 0].imshow(get_fields(field3_sum, norm_func), cmap='RdBu', interpolation='none',
    vmin=0.0, vmax=1.1,
    extent=(min_coords.x, max_coords.x, min_coords.y, max_coords.y), animated = True)
fig.colorbar(im3, ax=axes[1, 0])
axes[1, 0].set_xlabel("x")
axes[1, 0].set_ylabel("y")
axes[1, 0].set_title("field0 + field3")

# show field4
im4 = axes[1, 1].imshow(get_fields(field4_sum, norm_func), cmap='RdBu', interpolation='none',
    vmin=0.0, vmax=1.1,
    extent=(min_coords.x, max_coords.x, min_coords.y, max_coords.y), animated = True)
fig.colorbar(im2, ax=axes[1, 1])
axes[1, 1].set_xlabel("x")
axes[1, 1].set_ylabel("y")
axes[1, 1].set_title("field0 + field4")


def update_fig(*args):
    update_data()
    im1.set_array(get_fields(field1_sum, norm_func))
    im2.set_array(get_fields(field2_sum, norm_func))
    im3.set_array(get_fields(field3_sum, norm_func))
    im4.set_array(get_fields(field4_sum, norm_func))
    return im1, im2, im3, im4,
    
ani = animation.FuncAnimation(fig, update_fig, interval=10, blit=True)

fig.tight_layout()

plt.show()
