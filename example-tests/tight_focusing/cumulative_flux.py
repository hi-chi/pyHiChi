import sys
sys.path.append("../../bin")
import pyHiChi as hichi
import matplotlib.pyplot as plt
import numpy as np


def show(field, title, min_coords, max_coords):

    def get_fields(field, min_coords, max_coords, N):
        x = np.linspace(min_coords.x, max_coords.x, N)
        y = np.linspace(min_coords.y, max_coords.y, N)
        res = np.zeros(shape=(N,N))
        for ix in range(N):
            for iy in range(N):
                res[N - iy - 1, ix] = field.get_E(x[ix], y[iy], 0.0).norm()
        return res

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.imshow(get_fields(field, min_coords, max_coords, 256), cmap='RdBu', interpolation='none',
        extent=(min_coords.x, max_coords.x, min_coords.y, max_coords.y)
    )
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(title)
    plt.show()
    plt.close(fig)


wavelength = 1.0
R0 = 16*wavelength
pulselength = 2.0*wavelength
f_number = 0.3

tight_focusing_conf = hichi.TightFocusingField(f_number=f_number,
                                               R0=R0,
                                               wavelength=wavelength,
                                               pulselength=pulselength,
                                               total_power=hichi.c,
                                               edge_smoothing_angle=0.1
                                              )
FACTOR = 1.0
NX = int(FACTOR*200)
NY = int(FACTOR*200)
NZ = int(FACTOR*200)
grid_size = hichi.Vector3d(NX, NY, NZ)

focus_time_step = R0 / hichi.c                    

min_coords = hichi.Vector3d(-20*wavelength, -20*wavelength, -20*wavelength)
max_coords = hichi.Vector3d(20*wavelength, 20*wavelength, 20*wavelength)
grid_step = (max_coords - min_coords) / grid_size


# ------- initialize field -------
field = hichi.PSATDPoissonField(grid_size, min_coords, grid_step, focus_time_step)
field.set(tight_focusing_conf)
#show(field, "Start field", min_coords, max_coords)


# ------- get field in focus -------
field.update_fields()
#show(field, "Focused field", min_coords, max_coords)


# ------- zoom area near focus -------
focus_area_size = pulselength / 2

zoomed_min_coords = hichi.Vector3d(-focus_area_size, -focus_area_size, -focus_area_size)
zoomed_grid_size = grid_size*(zoomed_min_coords/min_coords)
zoomed_max_coords = zoomed_min_coords + zoomed_grid_size*grid_step
nx = int(zoomed_grid_size.x)
ny = int(zoomed_grid_size.y)
nz = int(zoomed_grid_size.z)

zoomed_field = field.zoom(
    min_coord=zoomed_min_coords,
    zoomed_grid_size=zoomed_grid_size,
    zoomed_grid_step=grid_step
)
#show(zoomed_field, "Zoomed field", zoomed_min_coords, zoomed_max_coords)


# ------- get f(t, y, z)=pulse(t, xi, y, z) -------

ntimes = 512
time_start = focus_time_step - focus_area_size / hichi.c
time_final = focus_time_step + focus_area_size / hichi.c
time_step = (time_final - time_start) / ntimes
times = np.arange(time_start, time_final, time_step)

xi_pos = -pulselength/4

#ex_yzt = np.zeros(shape=(ntimes, ny, nz), dtype='float64')
ey_yzt = np.zeros(shape=(ntimes, ny, nz), dtype='float64')
ez_yzt = np.zeros(shape=(ntimes, ny, nz), dtype='float64')
#bx_yzt = np.zeros(shape=(ntimes, ny, nz), dtype='float64')
by_yzt = np.zeros(shape=(ntimes, ny, nz), dtype='float64')
bz_yzt = np.zeros(shape=(ntimes, ny, nz), dtype='float64')

zoomed_field.advance(-focus_area_size / hichi.c)
zoomed_field.change_time_step(time_step)

for i in range(ntimes):
    args = {"x_pos":xi_pos,
        "y_min":zoomed_min_coords.y, "y_max":zoomed_max_coords.y, "y_size":int(zoomed_grid_size.y), 
        "z_min":zoomed_min_coords.z, "z_max":zoomed_max_coords.z, "z_size":int(zoomed_grid_size.z) 
    }
    #ex_yzt[i, :, :] = zoomed_field.get_Ez_yz_plane(**args)
    ey_yzt[i, :, :] = zoomed_field.get_Ey_yz_plane(**args)
    ez_yzt[i, :, :] = zoomed_field.get_Ez_yz_plane(**args)
    #bx_yzt[i, :, :] = zoomed_field.get_Bx_yz_plane(**args)
    by_yzt[i, :, :] = zoomed_field.get_By_yz_plane(**args)
    bz_yzt[i, :, :] = zoomed_field.get_Bz_yz_plane(**args)
    zoomed_field.update_fields()


# ------- compute integral -------

I = np.sum(ey_yzt*bz_yzt - ez_yzt*by_yzt, axis=0)
plt.matshow(I,
    extent=(zoomed_min_coords.y, zoomed_max_coords.y, zoomed_min_coords.z, zoomed_max_coords.z)
)
plt.title("I(xi, y, z)")
plt.xlabel("$y/\\lambda$")
plt.ylabel("$z/\\lambda$")
plt.show()

