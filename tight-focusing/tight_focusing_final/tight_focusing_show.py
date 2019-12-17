import sys
import os
sys.path.append("../../build/src/pyHiChi/Release")
import pyHiChi as hichi
import math as ma
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

Nx = 100
Ny = 100

minCoords = hichi.vector3d(0, 0, 0)
maxCoords = hichi.vector3d(1, 1, 1)

x = np.arange(minCoords.x, maxCoords.x, (maxCoords.x - minCoords.x)/Nx)
y = np.arange(minCoords.y, maxCoords.y, (maxCoords.y - minCoords.y)/Ny)

dirResult = "./pictures/"

fig = plt.figure()

vmin = 0
vmax = 0

def createDir(dir):
    if (os.path.exists(dir)): 
        for (dirpath, dirnames, filenames) in os.walk(dir):
            for file in filenames:
                os.remove(dir + file)
    else: os.mkdir(dir)


def createFieldAx(a, b, c, text):
    field = np.zeros(shape=(Ny,Nx))
    ax = fig.add_subplot(a, b, c)
    ax.title.set_text(text)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.tick_params(axis='both', which='major', labelsize=10)
    im = ax.imshow(field, cmap='RdBu', interpolation='none', extent=(minCoords.x, maxCoords.x, minCoords.y, maxCoords.y),\
        animated = True, aspect='auto', vmax=vmax, vmin=vmin)
    fig.colorbar(im, ax=ax)
    return im


def initVisual(minCoords_, maxCoords_, Nx_ = 100, Ny_ = 100, vmax_ = 5*10**8, vmin_ = 0):
   global x, y, Nx, Ny, minCoords, maxCoords, fig, vmax, vmin
   Nx = Nx_
   Ny = Ny_
   minCoords = minCoords_
   maxCoords = maxCoords_
   vmin = vmin_
   vmax = vmax_
   x = np.arange(minCoords.x, maxCoords.x, (maxCoords.x - minCoords.x)/Nx)
   y = np.arange(minCoords.y, maxCoords.y, (maxCoords.y - minCoords.y)/Ny)


def show(grid, update, maxIter=150):

    im = createFieldAx(1, 1, 1, "|E|")
    
    createDir(dirResult)

    def getFields():
        field = np.zeros(shape=(Ny,Nx))
        for ix in range(Nx):
            for iy in range(Ny):
                coordXY = hichi.vector3d(x[ix], y[iy], 0.0)
                E = grid.getE(coordXY)
                field[iy, ix] = E.norm()
        return field
    
    def animate(i):
        if (i > maxIter):
	        exit()
        print("\r %d" % (i), end = "")
        plt.savefig(dirResult + 'field%04d.png' % (i))
        update()
        field = getFields()
        im.set_array(field)
        return im,   

    ani = animation.FuncAnimation(fig, animate, interval=1, blit=True) 

    plt.show()



      
    