import sys
sys.path.append("../bin/")
import pyHiChi as hichi
import numpy as np
import math as ma

def valueEx(x, y, z):
    Ex = 0  #for x or y
    #Ex=np.sin(z) #for z
    return Ex

def valueEy(x, y, z):
    #Ey = 0 #for y or z
    #Ey = np.sin(x) #for x
    Ey = np.sin(x - z) #for xz
    return Ey

def valueEz(x, y, z):
    Ez = 0 #for x or z or xz
    #Ez = np.sin(y) #for y
    return Ez


def valueBx(x, y, z):
    #Bx = 0  #for x or z
    #Bx = np.sin(y) #for y
    Bx = np.sin(x - z)/np.sqrt(2) #for xz
    return Bx

def valueBy(x, y, z):
    By = 0  #for x or y or xz
    #By = np.sin(z) #for z
    return By

def valueBz(x, y, z):
    #Bz = 0  #for y or z
    #Bz = np.sin(x) #for x
    Bz = np.sin(x - z)/np.sqrt(2) #for xz
    return Bz


grid_size = hichi.Vector3d(20, 20, 20)
min_coords = hichi.Vector3d(0.0, 0.0, 0.0)
max_coords = hichi.Vector3d(2*np.pi, 2*np.pi, 2*np.pi)
grid_step = (max_coords - min_coords) / grid_size
time_step = 1e-14

field = hichi.YeeField(grid_size, min_coords, grid_step, time_step)
field.set_E(valueEx, valueEy, valueEz)
field.set_B(valueBx, valueBy, valueBz)

field.set_PML(0, 0, 0)
periodical_BC = hichi.PeriodicalBC(field)

#show
import matplotlib.pyplot as plt
import matplotlib.animation as animation

N = 50
eps = 0.0
x = np.arange(eps, 2*np.pi - eps, 2*(np.pi-eps)/N)
z = np.arange(eps, 2*np.pi - eps, 2*(np.pi-eps)/N)


def get_fields():
    global field, x, z, N
    #print(field)
    Ex = np.zeros(shape=(N, N))
    Ey = np.zeros(shape=(N, N))
    Ez = np.zeros(shape=(N, N))
    Bx = np.zeros(shape=(N, N))
    By = np.zeros(shape=(N, N))
    Bz = np.zeros(shape=(N, N))
    for ix in range(N):
        for iy in range(N):
            coord_xz = hichi.Vector3d(x[ix], 0.0, z[iy]) #for x or z or xz
            #coord_xz = hichi.Vector3d(x[ix], z[iy], 0.0) #for y or x
            
            E = field.get_E(coord_xz)
            Ex[ix, iy] = E.x
            Ey[ix, iy] = E.y
            Ez[ix, iy] = E.z
            B = field.get_B(coord_xz)
            Bx[ix, iy] = B.x
            By[ix, iy] = B.y
            Bz[ix, iy] = B.z
    return Ex, Ey, Ez, Bx, By, Bz

def update_data():
    for i in range(1000):
        field.update_fields()
    

(Ex, Ey, Ez, Bx, By, Bz) = get_fields()

fig, axes = plt.subplots(ncols=3, nrows=2)

im11 = axes[0, 0].imshow(Ex, cmap='RdBu', interpolation='none', extent=(0, 2*np.pi, 0, 2*np.pi), animated = True)
fig.colorbar(im11, ax=axes[0, 0])
axes[0, 0].set_title("Ex")
axes[0, 0].set_xlabel("x")
axes[0, 0].set_ylabel("z")

im12 = axes[0, 1].imshow(Ey, cmap='RdBu', interpolation='none', extent=(0, 2*np.pi, 0, 2*np.pi), animated = True)
fig.colorbar(im12, ax=axes[0, 1])
axes[0, 1].set_title("Ey")
axes[0, 1].set_xlabel("x")
axes[0, 1].set_ylabel("z")

im13 = axes[0, 2].imshow(Ez, cmap='RdBu', interpolation='none', extent=(0, 2*np.pi, 0, 2*np.pi), animated = True)
fig.colorbar(im13, ax=axes[0, 2])
axes[0, 2].set_title("Ez")
axes[0, 2].set_xlabel("x")
axes[0, 2].set_ylabel("z")

im21 = axes[1, 0].imshow(Bx, cmap='RdBu', interpolation='none', extent=(0, 2*np.pi, 0, 2*np.pi), animated = True)
fig.colorbar(im21, ax=axes[1, 0])
axes[1, 0].set_title("Bx")
axes[1, 0].set_xlabel("x")
axes[1, 0].set_ylabel("z")

im22 = axes[1, 1].imshow(By, cmap='RdBu', interpolation='none', extent=(0, 2*np.pi, 0, 2*np.pi), animated = True)
fig.colorbar(im22, ax=axes[1, 1])
axes[1, 1].set_title("By")
axes[1, 1].set_xlabel("x")
axes[1, 1].set_ylabel("z")

im23 = axes[1, 2].imshow(Bz, cmap='RdBu', interpolation='none', extent=(0, 2*np.pi, 0, 2*np.pi), animated = True)
fig.colorbar(im23, ax=axes[1, 2])
axes[1, 2].set_title("Bz")
axes[1, 2].set_xlabel("x")
axes[1, 2].set_ylabel("z")

i = 0

def update_fig(*args):
    global i
    update_data()
    (Ex, Ey, Ez, Bx, By, Bz) = get_fields()
    im11.set_array(Ex)
    im12.set_array(Ey)
    im13.set_array(Ez)
    im21.set_array(Bx)
    im22.set_array(By)
    im23.set_array(Bz)
    i = i + 1
    return im11, im12, im13, im21, im22, im23, 
    
ani = animation.FuncAnimation(fig, update_fig, interval=50, blit=True)

plt.tight_layout()

plt.show()


