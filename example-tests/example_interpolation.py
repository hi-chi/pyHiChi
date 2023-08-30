import sys
sys.path.append("../bin/")
import pyHiChi as hichi
import numpy as np

# ------------- model parameters -------------

N = 64
A = -1.0
B = 1.0
D = (B - A) / N
min_coords = hichi.Vector3d(-8*D, -8*D, A)
max_coords = hichi.Vector3d(8*D, 8*D, B)
grid_size = hichi.Vector3d(16, 16, N)

# field is a plane wave along z in 3d space

def value_Ex(x, y, z):
    Ex = np.sin(2*np.pi/(max_coords.z - min_coords.z)*z)
    return Ex

def value_Ey(x, y, z):
    Ey = 0
    return Ey

def value_Ez(x, y, z):
    Ez = 0
    return Ez

def value_Bx(x, y, z):
    Bx = 0
    return Bx

def value_By(x, y, z):
    By = np.sin(2*np.pi/(max_coords.z - min_coords.z)*z)
    return By

def value_Bz(x, y, z):
    Bz = 0
    return Bz

def analytical_solution(z, t):
    return np.sin(2*np.pi/(max_coords.z - min_coords.z)*(z - hichi.c*t))
    
grid_step = (max_coords - min_coords) / grid_size
time_step = grid_step.norm() / hichi.c * 0.1

# PSATD doesn't have numerical dispersion and shows very small numerical error for this task
field = hichi.PSATDField(grid_size, min_coords, grid_step, time_step)
field.set_E(value_Ex, value_Ey, value_Ez)
field.set_B(value_Bx, value_By, value_Bz)

def get_fields(a, b, n):
    # get 1d slice from numerical grid
    # interpolation is used if cells are not matched
    Ex = field.get_Ex_z_line(x_pos=0.0, y_pos=0.0,  # center coordinates
        z_min=a, z_max=b, z_size=n
    )
    return Ex
  
def get_analytical_fields(t, a, b, n):
    Z = np.linspace(a, b, n, endpoint=False)
    return analytical_solution(Z, t)

# ------------- show parameters -------------

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# If numerical grid size is great or equal to image construction grid size and multiply of it,
# then CIC, second order, fourth order interpolations show small error.
# However, it's just because of cell matching, not real interpolation error.

# image construction grid size = factor * numerical grid size
factors = [
    0.5,  # multiply, we do not see any interpolation error
    0.75, # not multiply, we can see interpolation error
    1.0,  # original size, we do not see any interpolation error
    2.0,  # greater, we can see interpolation error 
    4.0,  # greater, we can see interpolation error 
    8.0   # greater, we can see interpolation error 
]

interpolation_types = ["CIC", "TSC", "PCS", "second order", "fourth order"]
interpolation_linestyles = ["-", "--", "--", "--", "-."]
interpolation_funcs = [field.set_CIC_interpolation, field.set_TSC_interpolation, field.set_PCS_interpolation,
    field.set_second_order_interpolation, field.set_fourth_order_interpolation]

fig = plt.figure(layout="constrained")
gridspec = fig.add_gridspec(nrows=1, ncols=2, width_ratios=[1, 4])
gridspec_field = gridspec[0].subgridspec(2, 1)
gridspec_error = gridspec[1].subgridspec(3, 2)

title_fontsize = 10

N_EDGE_CELLS = 4  # temporary: interpolation near border may be incorrect, exclude edge cells
A_PLOT = A + N_EDGE_CELLS*grid_step.z
B_PLOT = B - N_EDGE_CELLS*grid_step.z
N_PLOT = N - 2*N_EDGE_CELLS

# create axes with EM field, numerical and analytical
def create_axis_fields(gs, labels, linestyles, title):
    ax = fig.add_subplot(gs)
    lines = []
    for label, linestyle in zip(labels, linestyles):
        lines.append(ax.plot([N_EDGE_CELLS, N_EDGE_CELLS + N_PLOT], [1, 1], linestyle, label=label)[0])
    ax.set_ylim((-1.0, 1.0))
    ax.set_title(title, fontsize=title_fontsize)
    ax.set_xlabel("x")
    ax.set_ylabel("field value")
    ax.grid()
    ax.legend()
    return lines

# create axes with error=abs(numerical - analytical)
def create_axis_errors(gs, labels, linestyles, title):
    ax = fig.add_subplot(gs)
    lines = []
    text = ""
    for label, linestyle in zip(labels, linestyles):
        lines.append(ax.plot([N_EDGE_CELLS, N_EDGE_CELLS + N_PLOT], [1, 1], linestyle, label=label)[0])
    ax.semilogy()
    ax.set_ylim((1e-14, 1e0))
    ax.set_title(title, fontsize=title_fontsize)
    ax.set_xlabel("numerical grid cell number")
    ax.set_ylabel("error")
    ax.grid()
    ax.legend(loc="upper left", ncol=len(labels))
    return lines

lines_field = []
lines_error = []

lines_field.extend(create_axis_fields(gridspec_field[0], ["numerical", "analytical"], ["-", "--"], "Field, grid size=%d" % (N)))
lines_field.extend(create_axis_errors(gridspec_field[1], ["|numerical-analytical|"], ["-"], "Error, grid size=%d" % (N)))
   
for row in range(3):
    for col in range(2):
        CUR_N = int(factors[row*2 + col] * N)
        lines_error.append(create_axis_errors(gridspec_error[row, col],
            interpolation_types, interpolation_linestyles, "Interpolation error, N=%d" % (CUR_N)))

# ------------- animation -------------

def update_fig(iteration):
    # field solver iteration
    field.update_fields()
    t = field.get_time()
    
    # show EM field, image construction grid and numerical grid are matched
    field.set_CIC_interpolation()
    
    X_TICKS = np.linspace(N_EDGE_CELLS, N - N_EDGE_CELLS, N - 2*N_EDGE_CELLS, endpoint=False)
    slice_non_edges = slice(N_EDGE_CELLS, -N_EDGE_CELLS)
    
    numerical_field = get_fields(min_coords.z, max_coords.z, N)
    lines_field[0].set_data(X_TICKS, numerical_field[slice_non_edges])
    
    analytical_field = get_analytical_fields(t, min_coords.z, max_coords.z, N)
    lines_field[1].set_data(X_TICKS, analytical_field[slice_non_edges])
    
    lines_field[2].set_data(X_TICKS, np.abs(numerical_field - analytical_field)[slice_non_edges])
        
    # show interpolation error for different image construction grid sizes  
    for j, func_set_interpolation, interpolation_type in zip(range(5), interpolation_funcs, interpolation_types):
        func_set_interpolation()
        
        for i, factor in enumerate(factors): 
            CUR_N_PLOT = int(factor*N_PLOT)  # current image construction grid size = factor * real grid size
            X_TICKS = np.linspace(N_EDGE_CELLS, N - N_EDGE_CELLS, CUR_N_PLOT, endpoint=False)
            
            interpolated_field = get_fields(A_PLOT, B_PLOT, CUR_N_PLOT)
            analytical_field = get_analytical_fields(t, A_PLOT, B_PLOT, CUR_N_PLOT)
            # show error depending on the current interpolation type
            lines_error[i][j].set_data(X_TICKS, np.abs(interpolated_field - analytical_field))
    
    return lines_field + [item for sublist in lines_error for item in sublist]
    
ani = animation.FuncAnimation(fig, update_fig, interval=50, blit=True)

plt.show()
