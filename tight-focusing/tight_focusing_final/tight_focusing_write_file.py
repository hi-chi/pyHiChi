import sys
import os
sys.path.append("../../build/src/pyHiChi/Release")
import pyHiChi as hichi
import numpy as np


def write(grid, update, minCoords, maxCoords, Nx = 300, Ny = 300, maxIter=160, dumpIter = 20, fileName = "res_x_%d.csv", dirResult = "./results/"):
    
    x = np.arange(minCoords.x, maxCoords.x, (maxCoords.x - minCoords.x)/Nx)
    y = np.arange(minCoords.y, maxCoords.y, (maxCoords.y - minCoords.y)/Ny)
    
    def getFields():
        field = np.zeros(shape=(Ny,Nx))
        for iy in range(Ny):
            for ix in range(Nx):
                coordXY = hichi.vector3d(x[ix], y[iy], 0.0)
                E = grid.getE(coordXY)
                field[iy, ix] = E.norm()
        return field
    
    for iter in range(maxIter + 1):
        print("\r %d" % iter),
        update()
        field = getFields()
        if (iter % dumpIter == 0):
            with open(dirResult + fileName % iter, "w") as file:
                for j in range(Ny):
                    for i in range(Nx):
                        file.write(str(field[j, i])+";")
                    file.write("\n")  
    print


def writeOX(grid, update, minCoords, maxCoords, Nx = 300, maxIter=160, dumpIter = 20, fileName = "res_x_%d.csv",\
    dirResult = "./results/", ifWriteZeroIter = True):
    
    x = np.arange(minCoords.x, maxCoords.x, (maxCoords.x - minCoords.x)/Nx)
    
    def getFields():
        field = np.zeros(shape=(Nx))
        for ix in range(Nx):
            coordXY = hichi.vector3d(x[ix], 0.0, 0.0)
            E = grid.getE(coordXY)
            field[ix] = E.norm()
        return field
    
    for iter in range(maxIter + 1):
        print("\r %d" % iter),
        update()
        field = getFields()
        if ((iter == 0 and ifWriteZeroIter) or (iter != 0 and iter % dumpIter == 0)):
            with open(dirResult + fileName % iter, "w") as file:
                for i in range(Nx):
                    file.write(str(field[i])+"\n")  
    print

