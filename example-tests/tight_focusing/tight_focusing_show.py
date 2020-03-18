import sys
import os
import pyHiChi as hichi
import math as ma
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import hichi_primitives

Nx = None
Ny = None

minCoords = None
maxCoords = None

x = None
y = None

vmin = None
vmax = None

fig = None
matplotlib.rcParams.update({"font.size" : 17})


def createFieldAx(a, b, c, text):
    field = np.zeros(shape=(Ny,Nx))
    ax = fig.add_subplot(a, b, c)
    ax.title.set_text(text)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.tick_params(axis='both', which='major')
    im = ax.imshow(field, cmap='RdBu', interpolation='none', extent=(minCoords.x, maxCoords.x, minCoords.y, maxCoords.y),\
        animated = True, aspect='auto', vmax=vmax, vmin=vmin)
    fig.colorbar(im, ax=ax)
    return im


# should be called first
def initVisual(minCoords_, maxCoords_, Nx_ = 300, Ny_ = 300, vmax_ = 0.5, vmin_ = 0):
   global x, y, Nx, Ny, minCoords, maxCoords, fig, vmax, vmin
   if (not fig == None): fig.clf()
   fig = plt.figure()
   Nx = Nx_
   Ny = Ny_
   minCoords = minCoords_
   maxCoords = maxCoords_
   vmin = vmin_
   vmax = vmax_
   x = np.arange(minCoords.x, maxCoords.x, (maxCoords.x - minCoords.x)/Nx)
   y = np.arange(minCoords.y, maxCoords.y, (maxCoords.y - minCoords.y)/Ny)


def getFields(grid):
    field = np.zeros(shape=(Ny,Nx))
    for ix in range(Nx):
        for iy in range(Ny):
            coordXY = hichi.vector3d(x[ix], y[iy], 0.0)
            E = grid.getE(coordXY)
            field[iy, ix] = E.norm()
    return field


# animate
def animate(grid, update, maxIter=160):

    im = createFieldAx(1, 1, 1, "|E|")
    fig.tight_layout()
    
    def animate_(i):
        if (i > maxIter):
            exit()
        update()
        field = getFields(grid)
        im.set_array(field)
        return im,   

    ani = animation.FuncAnimation(fig, animate_, interval=1, blit=True) 

    plt.show()


# save one picture
def savePicture(grid, dirResult = "./pictures/", name = "/res.png", dpi = 500):

    im = createFieldAx(1, 1, 1, "|E|")
    fig.tight_layout()
    
    field = getFields(grid)
    im.set_array(field)
    plt.savefig(dirResult + "/" + name, dpi=500)


# run and save pictures
def savePictures(grid, update, maxIter=160, dirResult = "./pictures/", dpi = 500):
    
    im = createFieldAx(1, 1, 1, "|E|")
    fig.tight_layout()
    
    hichi_primitives.createDir(dirResult)
    
    field = getFields(grid)
    im.set_array(field)
    plt.savefig(dirResult + '/field0000.png', dpi=500)

    for i in range(maxIter):
        update()
        field = getFields(grid)
        im.set_array(field)
        plt.savefig(dirResult + '/field%04d.png' % (i+1), dpi=dpi)

      
    