import sys
import os
import pyHiChi as hichi
import numpy as np


# write |E| to file, XOY plane
def writeXOY(grid, minCoords, maxCoords, Nx = 300, Ny = 300,
    fileName = "res_xoy.csv", dirResult = "./results/",
    x = None, y = None):
    
    if (type(x) == None): x = np.arange(minCoords.x, maxCoords.x, (maxCoords.x - minCoords.x)/Nx)
    if (type(y) == None): y = np.arange(minCoords.y, maxCoords.y, (maxCoords.y - minCoords.y)/Ny)
    
    def getFields():
        field = np.zeros(shape=(Ny,Nx))
        for iy in range(Ny):
            for ix in range(Nx):
                coordXY = hichi.vector3d(x[ix], y[iy], 0.0)
                E = grid.getE(coordXY)
                field[iy, ix] = E.norm()
        return field
        
    def write():
        field = getFields(grid)
        with open(dirResult + "/" + fileName, "w") as file:
            for j in range(Ny):
                for i in range(Nx):
                    file.write(str(field[j, i])+";")
                file.write("\n")
                
    write()


# run and write |E| to file, XOY plane
def runAndWriteXOY(grid, update, minCoords, maxCoords, Nx = 300, Ny = 300, maxIter=160, dumpIter = 20,
    fileName = "res_x_%d.csv", dirResult = "./results/", ifWriteZeroIter = True):
    
    if (ifWriteZeroIter):
        writeXOY(grid, minCoords, maxCoords, Nx, Ny, fileName % 0, dirResult, x, y)
    
    for iter in range(1, maxIter + 1):
        update()
        if (iter % dumpIter == 0):
            writeXOY(grid, minCoords, maxCoords, Nx, Ny, fileName % iter, dirResult, x, y)



# write |E| to file, OX axis
def writeOX(grid, minCoords, maxCoords, Nx = 300,
    fileName = "res_x_%d.csv", dirResult = "./results/",
    x = None):
    
    if (type(x) == None): x = np.arange(minCoords.x, maxCoords.x, (maxCoords.x - minCoords.x)/Nx)
    
    def getFields():
        field = np.zeros(shape=(Nx))
        for ix in range(Nx):
            coordXY = hichi.vector3d(x[ix], 0.0, 0.0)
            E = grid.getE(coordXY)
            field[ix] = E.norm()
        return field
        
    def write():
        field = getFields()
        with open(dirResult + "/" + fileName, "w") as file:
            for i in range(Nx):
                file.write(str(field[i])+"\n")
        
    write()


# run and write |E| to file, OX axis
def runAndWriteOX(grid, update, minCoords, maxCoords, Nx = 300, maxIter=160, dumpIter = 20, fileName = "res_x_%d.csv",
    dirResult = "./results/", ifWriteZeroIter = True):
    
    x = np.arange(minCoords.x, maxCoords.x, (maxCoords.x - minCoords.x)/Nx)
        
    if (ifWriteZeroIter):
        writeOX(grid, minCoords, maxCoords, Nx, fileName % 0, dirResult, x = x)

    for iter in range(1, maxIter + 1):
        update()
        if (iter % dumpIter == 0):
            writeOX(grid, minCoords, maxCoords, Nx, fileName % iter, dirResult, x = x)

