import sys
import os
import pyHiChi as hichi
import numpy as np


# run and write |E| to file, XOY plane
def writeXOY(grid, update, minCoords, maxCoords, Nx = 300, Ny = 300, maxIter=160, dumpIter = 20,
    fileName = "res_x_%d.csv", dirResult = "./results/", ifWriteZeroIter = True):
    
    x = np.arange(minCoords.x, maxCoords.x, (maxCoords.x - minCoords.x)/Nx)
    y = np.arange(minCoords.y, maxCoords.y, (maxCoords.y - minCoords.y)/Ny)
    
    def getFields(grid):
        field = np.zeros(shape=(Ny,Nx))
        for iy in range(Ny):
            for ix in range(Nx):
                coordXY = hichi.vector3d(x[ix], y[iy], 0.0)
                E = grid.getE(coordXY)
                field[iy, ix] = E.norm()
        return field
        
    def write(grid, iter):
        field = getFields(grid)
        with open(dirResult + "/" + fileName % iter, "w") as file:
            for j in range(Ny):
                for i in range(Nx):
                    file.write(str(field[j, i])+";")
                file.write("\n")
                
    if (ifWriteZeroIter):
        write(grid, 0)
    
    for iter in range(1, maxIter + 1):
        update()
        if (iter % dumpIter == 0):
            write(grid, iter)


# run and write |E| to file, OX axis
def writeOX(grid, update, minCoords, maxCoords, Nx = 300, maxIter=160, dumpIter = 20, fileName = "res_x_%d.csv",\
    dirResult = "./results/", ifWriteZeroIter = True):
    
    x = np.arange(minCoords.x, maxCoords.x, (maxCoords.x - minCoords.x)/Nx)
    
    def getFields(grid):
        field = np.zeros(shape=(Nx))
        for ix in range(Nx):
            coordXY = hichi.vector3d(x[ix], 0.0, 0.0)
            E = grid.getE(coordXY)
            field[ix] = E.norm()
        return field
        
    def write(grid, iter):
        field = getFields(grid)
        with open(dirResult + "/" + fileName % (iter), "w") as file:
            for i in range(Nx):
                file.write(str(field[i])+"\n")
        
    if (ifWriteZeroIter):
        write(grid, 0)

    for iter in range(1, maxIter + 1):
        update()
        if (iter % dumpIter == 0):
            write(grid, iter)

