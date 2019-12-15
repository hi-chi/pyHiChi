'''
import sys
sys.path.append("C:/Codes/hichi/pyHiChi-master/build/visual_studio/src/pyHiChi/Release")
sys.path.append("C:/Codes/hichi/draft_test1")
from math import *
import numpy as np
import os
import pyHiChi as hichi
from hichi_primitives import *
import matplotlib.pyplot as plt

from matplotlib.colors import ListedColormap, LinearSegmentedColormap

print("show")

def getFields(x, y, Ex, Ey, E_abs, field_):
    E = hichi.vector3d(0, 0, 0)
    B = hichi.vector3d(0, 0, 0)
    E_max = 0
    Ny = y.shape[0]
    for ix in range(x.shape[0]):
        for iy in range(y.shape[0]):
            coordXY = hichi.vector3d(x[ix], y[iy], 0.0)
            field_.get(coordXY, E, B)
            Ex[Ny - 1 - iy, ix] = E.x
            Ey[Ny - 1 - iy, ix] = E.y
            E_abs[Ny - 1 - iy, ix] = E.norm()
            if(E_max < E.norm()): E_max = E.norm()
            #print(x[ix], y[iy], ' ->', E.x, E.y, E.z, '\n')
    return E_max



def showGrid_(grid, minX, maxX, minY, maxY, field_, showGridCond):
    grid[:][:] = 0
    it = field_.getIterator()
    it.selectPositions(showGridCond)
    it.begin()
    rr = hichi.vector3d(0, 0, 0)
    while (it.next(rr)):
        x = round(grid.shape[0] * (rr.x - minX) / (maxX - minX))
        y = round(grid.shape[1] * (rr.y - minY) / (maxY - minY))
        if(x < grid.shape[0] and x >= 0 and y < grid.shape[1] and y >= 0): grid[x][y] = 1


def field(field_, minX, maxX, minY, maxY, Nx, Ny, it = 0, interval = 8, showGridCond = nullFunc3):
    x = np.arange(minX, maxX, (maxX - minX) / Nx)
    y = np.arange(minY, maxY, (maxY - minY) / Ny)

    X = np.zeros(shape=(Ny, Nx))  # variable on the plot
    Y = np.zeros(shape=(Ny, Nx))  # variable on the plot
    Ex = np.zeros(shape=(Ny, Nx))
    Ey = np.zeros(shape=(Ny, Nx))
    E_abs = np.zeros(shape=(Ny, Nx))
    grid = np.zeros(shape=(Ny, Nx))

    for ix in range(x.shape[0]):
        for iy in range(y.shape[0]):
            X[Ny - 1 - iy, ix] = x[ix] / 1e-4
            Y[Ny - 1 - iy, ix] = y[iy] / 1e-4

    E_max = getFields(x, y, Ex, Ey, E_abs, field_)
    Ex = Ex/E_max
    Ey = Ey/E_max
    fig, ax = plt.subplots()
    ax.set_title('Electric field $\mathbf{E}$,  normalized to $E_{max}(t = 0) =$ ' + '{:.2e}'.format(E_max))
    ax.set_aspect(1.0)
    ax.set_xlabel('x [$\mu$m]')
    ax.set_ylabel('y [$\mu$m]')
    #im = ax.imshow(E_abs, cmap='YlOrBr', interpolation='none', extent=(minX/1e-4, maxX/1e-4, minY/1e-4, maxY/1e-4), vmax=E_max, vmin=0)
    showGrid_(grid, minX, maxX, minY, maxY, field_, showGridCond)
    im = ax.imshow(grid, cmap='YlOrBr', interpolation='none', extent=(minX / 1e-4, maxX / 1e-4, minY / 1e-4, maxY / 1e-4), vmax=2, vmin=0)
    i = interval
    col = plt.get_cmap('YlOrBr', 256)
    newcmp = ListedColormap(col(np.linspace(0, 1.5, 256)))
    q = ax.quiver(X[::i, ::i], Y[::i, ::i], Ex[::i, ::i], Ey[::i, ::i], E_abs[::i, ::i], pivot='mid', cmap=newcmp, scale=x.shape[0] / 8)

    #list_ = list()
    #radius = 0.001
    #if(showGrid):
    #    showGrid_(list_, ax, field_, radius, showGridCond)
    #for item in list_:
    #    item.remove()
    plt.savefig('field0.png', dpi=200)
    if(it > 1):
        dir = os.getcwd() + '/field'
        if not os.path.exists(dir):
            os.makedirs(dir)

        for j in np.arange(it):
            getFields(x, y, Ex, Ey, E_abs, field_)
            im.set_array(E_abs)
            Ex = Ex / E_max
            Ey = Ey / E_max
            q.set_UVC(Ex[::i, ::i], Ey[::i, ::i], E_abs[::i, ::i])
            list_ = list()
            #if (showGrid):
            #    showGrid_(list_, ax, field_, radius, showGridCond)
            plt.savefig((dir + '/field' + format(j, '04d') + '.png'), dpi=200)
            #for item in list_:
            #    item.remove()
            print('time = ', field_.getTime(), '\n')
            field_.advance()
'''
import sys
sys.path.append("C:/Users/Elena/Documents/Visual Studio 2017/Projects/plasma/hi-chi/pyHiChi/build/src/pyHiChi/Release")
#sys.path.append("C:/Codes/hichi/draft_test1")
from math import *
import numpy as np
import os
import pyHiChi as hichi
from hichi_primitives import *
import matplotlib.pyplot as plt

from matplotlib.colors import ListedColormap, LinearSegmentedColormap

print("show")

def getFields(x, y, Ex, Ey, E_abs, field_):
    E = hichi.vector3d(0, 0, 0)
    B = hichi.vector3d(0, 0, 0)
    E_max = 0
    Ny = y.shape[0]
    for ix in range(x.shape[0]):
        for iy in range(y.shape[0]):
            coordXY = hichi.vector3d(x[ix], y[iy], 0.0)
            field_.get(coordXY, E, B)
            Ex[Ny - 1 - iy, ix] = E.x
            Ey[Ny - 1 - iy, ix] = E.y
            E_abs[Ny - 1 - iy, ix] = E.norm()
            if(E_max < E.norm()): E_max = E.norm()
            #print(x[ix], y[iy], ' ->', E.x, E.y, E.z, '\n')
    return E_max

def showGrid_(list_, ax, field_, r, showGridCond):
    for item in list_:
        item.remove()
    del list_[:]
    it = field_.getIterator()
    it.selectPositions(showGridCond)
    it.begin()
    rr = hichi.vector3d(0, 0, 0)
    while (it.next(rr)):
        #print(rr.x, rr.y, '\n')
        circle = plt.Circle((rr.x / 1e-4, rr.y / 1e-4), r, color='k')
        list_.append(circle)
        ax.add_artist(circle)

def showGrid_1(grid, minX, maxX, minY, maxY, field_, showGridCond):
    grid[:][:] = 0
    it = field_.getIterator()
    it.selectPositions(showGridCond)
    it.begin()
    rr = hichi.vector3d(0, 0, 0)
    while (it.next(rr)):
        x = round(grid.shape[0] * (rr.x - minX) / (maxX - minX))
        y = round(grid.shape[1] * (rr.y - minY) / (maxY - minY))
        if(x < grid.shape[0] and x >= 0 and y < grid.shape[1] and y >= 0): grid[y][x] = 1



def field(field_, minX, maxX, minY, maxY, Nx, Ny, it = 0, interval = 8, showGridCond = nullFunc3):
    x = np.arange(minX, maxX, (maxX - minX) / Nx)
    y = np.arange(minY, maxY, (maxY - minY) / Ny)

    X = np.zeros(shape=(Ny, Nx))  # variable on the plot
    Y = np.zeros(shape=(Ny, Nx))  # variable on the plot
    Ex = np.zeros(shape=(Ny, Nx))
    Ey = np.zeros(shape=(Ny, Nx))
    E_abs = np.zeros(shape=(Ny, Nx))
    grid = np.zeros(shape=(round(Ny/1), round(Nx/1)))

    for ix in range(x.shape[0]):
        for iy in range(y.shape[0]):
            X[Ny - 1 - iy, ix] = x[ix] / 1e-4
            Y[Ny - 1 - iy, ix] = y[iy] / 1e-4

    E_max = getFields(x, y, Ex, Ey, E_abs, field_)
    Ex = Ex/E_max
    Ey = Ey/E_max
    fig, ax = plt.subplots()
    E_max = E_max*10
    #ax.set_title('Electric field $\mathbf{E}$,  normalized to ' + '{:.2e}'.format(E_max))
    ax.set_title('Electric field strength')
    ax.set_aspect(1.0)
    ax.set_xlabel('x [$\mu$m]')
    ax.set_ylabel('y [$\mu$m]')
    im = ax.imshow(E_abs, cmap='YlOrBr', interpolation='none', extent=(minX/1e-4, maxX/1e-4, minY/1e-4, maxY/1e-4), vmax=E_max, vmin=0)
    showGrid_1(grid, minX, maxX, minY, maxY, field_, showGridCond)
    xg = np.arange(minX/1e-4, maxX/1e-4, ((maxX - minX)/1e-4) / grid.shape[0])
    yg = np.arange(minY/1e-4, maxY/1e-4, ((maxY - minY)/1e-4) / grid.shape[1])
    Xg, Yg = np.meshgrid(xg, yg)
    img = ax.contour(Xg, Yg, grid, 1, linewidths=0.25, colors='k')
    i = interval
    col = plt.get_cmap('YlOrBr', 256)
    newcmp = ListedColormap(col(np.linspace(0, 1.5, 256)))
    #q = ax.quiver(X[::i, ::i], Y[::i, ::i], Ex[::i, ::i], Ey[::i, ::i], E_abs[::i, ::i], pivot='mid', cmap=newcmp, scale=x.shape[0] / 8)
    list_ = list()
    radius = 0.001
    #showGrid_(list_, ax, field_, radius, showGridCond)

    plt.savefig('field0.png', dpi=300)
    if(it > 1):
        dir = os.getcwd() + '/field'
        if not os.path.exists(dir):
            os.makedirs(dir)

        for j in np.arange(it):
            getFields(x, y, Ex, Ey, E_abs, field_)
            im.set_array(E_abs)
            showGrid_1(grid, minX, maxX, minY, maxY, field_, showGridCond)
            for coll in img.collections:
                plt.gca().collections.remove(coll)
            img = ax.contour(Xg, Yg, grid, 1, linewidths=0.25, colors='k')
            Ex = Ex / E_max
            Ey = Ey / E_max
            #q.set_UVC(Ex[::i, ::i], Ey[::i, ::i], E_abs[::i, ::i])
            #showGrid_(list_, ax, field_, radius, showGridCond)
            plt.savefig((dir + '/field' + format(j, '04d') + '.png'), dpi=300)
            print('time = ', field_.getTime(), '\n')
            field_.advance()










