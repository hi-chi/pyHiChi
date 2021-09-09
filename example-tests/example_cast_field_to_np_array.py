import sys
sys.path.append("../bin")
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
    return np.exp(-x**2-y**2)*np.sin(3*x)


field = hichi.PSATDField(grid_size, min_coords, grid_step, time_step)
field.set_E(null_value, field_value, null_value)
field.set_B(null_value, null_value, field_value)
scalar_field = field.get_Ey()

# makes 3d np.array of Ey field component
scalar_field_arr = np.array(scalar_field)
print(scalar_field_arr)
print(scalar_field_arr.shape)  # shape of scalar field can be not equal to grid_size
                               # there may be more elements in the last dimension
                               
scalar_field_arr[N//2, N//2, 0] = 10000000.0
print(field.get_E(0.0, 0.0, 0.0).y)  # demonstrates that scalar_field_arr is the copy of real field

# makes 3d np.array of Ey field component without copying (be careful!)
scalar_field_arr = np.array(scalar_field, copy=False)
scalar_field_arr[N//2, N//2, 0] = 10000000.0  # rewrites value in field
print(field.get_E(0.0, 0.0, 0.0).y)

