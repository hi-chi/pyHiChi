import sys
sys.path.append("../bin/")
import pyHiChi as hichi
import numpy as np


def field_value(x, y, z):
    return np.exp(-x**2-y**2-z**2)*np.cos(5*x)
    
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

N = 128

x_line = field.get_Ey_x_line(y_pos=0.0, z_pos=0.0, x_min=min_coords.x, x_max=max_coords.x, x_size=N)
y_line = field.get_Ey_y_line(x_pos=0.0, z_pos=0.0, y_min=min_coords.y, y_max=max_coords.y, y_size=N)
z_line = field.get_Ey_z_line(x_pos=0.0, y_pos=0.0, z_min=min_coords.z, z_max=max_coords.z, z_size=N)
print(type(x_line))
print(x_line.shape)

xy_plane = field.get_Ey_xy_plane(z_pos=0.0, x_min=min_coords.x, x_max=max_coords.x, x_size=N,
    y_min=min_coords.y, y_max=max_coords.y, y_size=N)
xz_plane = field.get_Ey_xz_plane(y_pos=0.0, x_min=min_coords.x, x_max=max_coords.x, x_size=N,
    z_min=min_coords.z, z_max=max_coords.z, z_size=N)
yz_plane = field.get_Ey_yz_plane(x_pos=0.0, y_min=min_coords.y, y_max=max_coords.y, y_size=N,
    z_min=min_coords.z, z_max=max_coords.z, z_size=N)
print(type(xy_plane))
print(xy_plane.shape)

area_3d = field.get_Ey(x_min=min_coords.x, x_max=max_coords.x, x_size=N,
    y_min=min_coords.y, y_max=max_coords.y, y_size=N, z_min=min_coords.z, z_max=max_coords.z, z_size=N)
print(type(area_3d))
print(area_3d.shape)


# --------- show -------------


import matplotlib.pyplot as plt
import matplotlib.animation as animation

fig, axes = plt.subplots(ncols=3, nrows=2)


axes[0, 0].plot(np.linspace(min_coords.x, max_coords.x, N), x_line)
axes[0, 0].set_xlabel("x")
axes[0, 0].set_ylabel("Ey")
axes[0, 0].set_title("X line (y=0, z=0)")

axes[0, 1].plot(np.linspace(min_coords.y, max_coords.y, N), y_line)
axes[0, 1].set_xlabel("y")
axes[0, 1].set_ylabel("Ey")
axes[0, 1].set_title("Y line (x=0, z=0)")

axes[0, 2].plot(np.linspace(min_coords.z, max_coords.z, N), z_line)
axes[0, 2].set_xlabel("z")
axes[0, 2].set_ylabel("Ey")
axes[0, 2].set_title("Z line (x=0, y=0)")


im10 = axes[1, 0].imshow(xy_plane, cmap='RdBu', interpolation='none',
    extent=(min_coords.x, max_coords.x, min_coords.y, max_coords.y), animated = True)
fig.colorbar(im10, ax=axes[1, 0])
axes[1, 0].set_xlabel("x")
axes[1, 0].set_ylabel("y")
axes[1, 0].set_title("XY plane (z=0)")

im11 = axes[1, 1].imshow(xz_plane, cmap='RdBu', interpolation='none',
    extent=(min_coords.x, max_coords.x, min_coords.z, max_coords.z), animated = True)
fig.colorbar(im11, ax=axes[1, 1])
axes[1, 1].set_xlabel("x")
axes[1, 1].set_ylabel("z")
axes[1, 1].set_title("XZ plane (y=0)")

im12 = axes[1, 2].imshow(yz_plane, cmap='RdBu', interpolation='none',
    extent=(min_coords.y, max_coords.y, min_coords.z, max_coords.z), animated = True)
fig.colorbar(im12, ax=axes[1, 2])
axes[1, 2].set_xlabel("y")
axes[1, 2].set_ylabel("z")
axes[1, 2].set_title("YZ plane (x=0)")


fig.tight_layout()

plt.show()
