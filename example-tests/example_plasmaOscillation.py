import sys
sys.path.append("../bin")
import pyHiChi as hichi

import numpy as np
import math as ma
import numba
from numba import cfunc

from matplotlib import pyplot as plt

nx = 64
ny = 8
nz = 8
L = 1.0
dx = L / nx
dy = dx
dz = dx
Nip = 256
Np = 14
N = int(Nip * Np)
T0 = 0.01 * hichi.ELECTRON_MASS * hichi.c * hichi.c
D = T0 / (8 * hichi.pi * hichi.ELECTRON_CHARGE * hichi.ELECTRON_CHARGE * 0.25 * dx * dx)
wp = ma.sqrt(4.0 * hichi.pi * hichi.ELECTRON_CHARGE * hichi.ELECTRON_CHARGE * D / hichi.ELECTRON_MASS)
dt = 2 * (hichi.pi / wp) / Nip
Nc = 30
w = D * dx * dx * dx / Nc
A = 0.05

field_size = hichi.Vector3d(nx, ny, nz)
min_coords = hichi.Vector3d(0.0, 0.0, 0.0)
field_step = hichi.Vector3d(dx, dy, dz)

field = hichi.YeeField(field_size, min_coords, field_step, dt)

x = np.arange(min_coords.x, L, dx)

def valueEx(x, y, z):
    v1 = hichi.Vector3d(x, y, z)
    v2 = hichi.Vector3d(min_coords.x + v1.x, min_coords.y + v1.y, min_coords.z + v1.z)
    Ex = -2 * L * D * A * hichi.ELECTRON_CHARGE * ma.cos(2 * hichi.pi * v2.x / L)
    return Ex

def valueEy(x, y, z):
    Ey = 0
    return Ey

def valueEz(x, y, z):
    Ez = 0 
    return Ez

def valueBx(x, y, z):
    Bx = 0
    return Bx

def valueBy(x, y, z):
    By = 0  
    return By

def valueBz(x, y, z):
    Bz = 0
    return Bz

@cfunc("double(double,double,double)")
def DensityFunc(x, y, z):
    T0 = 0.01 * hichi.ELECTRON_MASS * hichi.c * hichi.c
    L = 1.0
    dx = L / 64.0
    D = T0 / (8 * hichi.pi * hichi.ELECTRON_CHARGE * hichi.ELECTRON_CHARGE * 0.25 * dx * dx)
    return D * (1 + 0.05 * ma.sin(2 * hichi.pi * x / L))


@cfunc("double(double,double,double)")
def InitialTemperatureFunc(x, y, z):
    T0 = 0.01 * hichi.ELECTRON_MASS * hichi.c * hichi.c
    return T0

def get_fields():
    global field, x, nx
    Ex = np.zeros(shape=(nx))
    for ix in range(nx):
        coord_x = hichi.Vector3d(x[ix], 0.0, 0.0)
        E = field.get_E(coord_x)
        Ex[ix] = E.x
    return Ex

def get_particle_density(particleArray):
    res = []
    for k in range(nx):
        electronCount = 0
        minCellCoords = min_coords.x + k * dx
        maxCellCoords = min_coords.x + (k + 1) * dx

        for j in range(particleArray.size()):
            if (particleArray[j].get_position().x >= minCellCoords) \
                and (particleArray[j].get_position().x <= maxCellCoords):
                electronCount = electronCount + 1
        electronDensity = electronCount * w  # density through plane
        res.append(electronDensity)
    return res


field.set_E(valueEx, valueEy, valueEz)
field.set_B(valueBx, valueBy, valueBz)
field.set_periodical_BC()

pusher = hichi.BorisPusher()
currentDeposition = hichi.DepositionCIC(field, dt)

J_BC = hichi.periodical_J_BC(field)
particle_BC = hichi.periodical_particle_BC(field)

p_array = hichi.ParticleArray()

particle_generator = hichi.ParticleGenerator(field)
particle_generator(p_array, DensityFunc.address, InitialTemperatureFunc.address, numba.f4(w), hichi.ELECTRON)
interpolation = hichi.InterpolationCIC(field)

fig = plt.figure()
def plotEx(index, title):
    ax = fig.add_subplot(1,2,index)
    ax.plot(x, get_fields())
    ax.set_ylabel("$E_x$")
    ax.set_xlim((min_coords.x, L))
    ax.set_xlabel("$x$")
    ax.grid()
    
    ax.set_title(title)

def plotParticleDens(density, index, title):
    ax = fig.add_subplot(1,2,index)
    ax.plot(x, density)
    ax.set_ylabel("$Particle Density$")
    ax.set_xlim((min_coords.x, L))
    ax.set_xlabel("$x$")
    ax.grid()
    
    ax.set_title(title)

nShow = 0
for i in range (N + 1):
    (Ex) = get_fields()
    # fdtd
    field.update_fields()

    # interpolate
    fields = []
    for j in range(p_array.size()):
        coords = hichi.Vector3d(p_array[j].get_position().x, p_array[j].get_position().y, p_array[j].get_position().z)
        E_x = interpolation.getExCIC(coords)
        E_y = interpolation.getEyCIC(coords)
        E_z = interpolation.getEzCIC(coords)
        E = hichi.Vector3d(E_x, E_y, E_z)
        B_x = interpolation.getBxCIC(coords)
        B_y = interpolation.getByCIC(coords)
        B_z = interpolation.getBzCIC(coords)
        B = hichi.Vector3d(B_x, B_y, B_z)
        fields.append(hichi.FieldValue(E, B))

    # pusher
    pusher(p_array, fields, dt)
    # periodical particle BC
    particle_BC.update(p_array)

    if i == nShow:
        plotEx(1, "E_x")
        plotParticleDens(get_particle_density(p_array), 2, "Particle Density")
    #current deposition
    currentDeposition(p_array)
    J_BC.update()

    if i % 16 == 0:
        print(i)

plt.tight_layout()
plt.show()
