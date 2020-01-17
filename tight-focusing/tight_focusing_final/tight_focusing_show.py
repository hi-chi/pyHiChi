import sys
import os
sys.path.append("../../build/src/pyHiChi/Release")
import pyHiChi as hichi
import math as ma
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.ticker as ticker
import numpy as np
import hichi_primitives

Nx = 300
Ny = 300

minCoords = hichi.vector3d(0, 0, 0)
maxCoords = hichi.vector3d(1, 1, 1)

x = np.arange(minCoords.x, maxCoords.x, (maxCoords.x - minCoords.x)/Nx)
y = np.arange(minCoords.y, maxCoords.y, (maxCoords.y - minCoords.y)/Ny)

fig = plt.figure()
matplotlib.rcParams.update({"font.size" : 17})

vmin = 0
vmax = 0


def createFieldAx(a, b, c, text):
    field = np.zeros(shape=(Ny,Nx))
    ax = fig.add_subplot(a, b, c)
    ax.title.set_text(text)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.tick_params(axis='both', which='major')
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.001))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.0005))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.001))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.0005))
    im = ax.imshow(field, cmap='RdBu', interpolation='none', extent=(minCoords.x, maxCoords.x, minCoords.y, maxCoords.y),\
        animated = True, aspect='auto', vmax=vmax, vmin=vmin)
    fig.colorbar(im, ax=ax)
    return im


def initVisual(minCoords_, maxCoords_, Nx_ = 300, Ny_ = 300, vmax_ = 5*10**8, vmin_ = 0):
   global x, y, Nx, Ny, minCoords, maxCoords, fig, vmax, vmin
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


def savePictures(grid, update, maxIter=160, dirResult = "./pictures/"):
    
    im = createFieldAx(1, 1, 1, "|E|")
    fig.tight_layout()
    
    hichi_primitives.createDir(dirResult)
    
    field = getFields(grid)
    im.set_array(field)
    plt.savefig(dirResult + 'field0000.png', dpi=1000)

    for i in range(maxIter):
        update()
        field = getFields(grid)
        im.set_array(field)
        plt.savefig(dirResult + 'field%04d.png' % (i+1), dpi=1000)

      
    