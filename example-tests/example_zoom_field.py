import sys
sys.path.append("../bin/")
import pyHiChi as hichi
import numpy as np

N = 128
min_coords = hichi.Vector3d(-10, -10, 0.0)
max_coords = hichi.Vector3d(10, 10, 0.0)
grid_size = hichi.Vector3d(N, N, 1)
d = (max_coords.x - min_coords.x) / grid_size.x
grid_step = hichi.Vector3d(d, d, d)
time_step = 2.5/hichi.c


def null_value(x, y, z):
    return 0.0

def field_value(x, y, z):
    return np.exp(-x**2-y**2)*np.sin(3*x)  # omega=3


field = hichi.PSATDField(grid_size, min_coords, grid_step, time_step)
field.set_E(null_value, field_value, null_value)
field.set_B(null_value, null_value, field_value)

# here we create zoomed copy of the field
zoomed_grid_size = hichi.Vector3d(N // 3, N // 3, 1)
zoomed_min_coords = (max_coords - min_coords) / 3 + min_coords
zoomed_max_coords = zoomed_min_coords + grid_step*zoomed_grid_size
zoomed_field = field.zoom(
    min_coord=zoomed_min_coords,
    zoomed_grid_size=zoomed_grid_size,
    zoomed_grid_step=grid_step
)

zoomed_field.update_fields()


# --------- show -------------

import matplotlib.pyplot as plt
import matplotlib.animation as animation

def get_fields(field, min_coords, max_coords, N):
    x = np.linspace(min_coords.x, max_coords.x, N)
    y = np.linspace(min_coords.y, max_coords.y, N)
    res = np.zeros(shape=(N,N))
    for ix in range(N):
        for iy in range(N):
            coord_xy = hichi.Vector3d(x[ix], y[iy], 0.0)
            res[N - iy - 1, ix] = field.get_E(coord_xy).norm()
    return res

fig = plt.figure()
ax1 = fig.add_subplot(1,2,1)
ax2 = fig.add_subplot(1,2,2)

image1 = get_fields(field, min_coords, max_coords, N)
im1 = ax1.imshow(image1, cmap='RdBu', interpolation='none',
    extent=(min_coords.x, max_coords.x, min_coords.y, max_coords.y)
)
fig.colorbar(im1, ax=ax1)
ax1.set_xlabel("x")
ax1.set_ylabel("y")
ax1.set_title("Start field")

image2 = get_fields(zoomed_field, zoomed_min_coords, zoomed_max_coords, N // 3)
im2 = ax2.imshow(image2, cmap='RdBu', interpolation='none',
    extent=(zoomed_min_coords.x, zoomed_max_coords.x, zoomed_min_coords.y, zoomed_max_coords.y)
)
fig.colorbar(im2, ax=ax2)
ax2.set_xlabel("x")
ax2.set_ylabel("y")
ax2.set_title("Zoomed advanced field")

fig.tight_layout()

plt.show()

