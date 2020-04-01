import pyHiChi as hichi
import numpy as np
import os
from hichi_primitives import Axis, Plane, Field, getCoordValue


class Writer:
    
    def __init__(self, grid, minCoords, maxCoords, dir="./"):
        self.grid = grid
        self.minCoords = minCoords
        self.maxCoords = maxCoords
        self.dir = os.path.join(os.getcwd(), dir)
        
    
    def saveFileInPlane(self, shape, plane=Plane.XOY, lastCoordinateValue=0.0,
                        field=Field.E, fieldCoord=Axis.X, norm=False,
                        nameFile="field.csv"):
        minCoords = (getCoordValue(self.minCoords, plane.value[0]), getCoordValue(self.minCoords, plane.value[1]))
        maxCoords = (getCoordValue(self.maxCoords, plane.value[0]), getCoordValue(self.maxCoords, plane.value[1]))
        coords0 = np.linspace(minCoords[0], maxCoords[0], shape[0])
        coords1 = np.linspace(minCoords[1], maxCoords[1], shape[1])         
        fields = self.getFieldPlane_((coords0, coords1), shape, plane, lastCoordinateValue,
            self.generateFuncGet_(field, fieldCoord, norm))
            
        with open(os.path.join(self.dir, nameFile), "w") as file:
            for iy in range(shape[1]):
                for ix in range(shape[0]):
                    file.write("%f;" % fields[iy, ix])
                file.write("\n")
            
        
    def saveFileInAxis(self, nPoints, axis=Axis.X, lastCoordinateValue=(0.0, 0.0),
                       field=Field.E, fieldCoord=Axis.X, norm=False,
                       nameFile="field.csv"):
        minCoords = getCoordValue(self.minCoords, axis)
        maxCoords = getCoordValue(self.maxCoords, axis)
        coords = np.linspace(minCoords, maxCoords, nPoints)
        fields = self.getFieldAxis_(coords, nPoints, axis, lastCoordinateValue,
            self.generateFuncGet_(field, fieldCoord, norm))
            
        with open(os.path.join(self.dir, nameFile), "w") as file:
            for ix in range(nPoints):
                file.write("%f\n" % fields[ix])
                
    
    def generateFuncGet_(self, field, fieldCoord, norm):
        func = None
        
        if field == Field.E:
            if norm:
                func = lambda self, coords: self.grid.getE(coords).norm()
            else:
                if fieldCoord == Axis.X:
                    func = lambda self, coords: self.grid.getE(coords).x
                elif fieldCoord == Axis.Y:
                    func = lambda self, coords: self.grid.getE(coords).y
                elif fieldCoord == Axis.Z:
                    func = lambda self, coords: self.grid.getE(coords).z
                    
        elif field == Field.B:
            if norm:
                func = lambda self, coords: self.grid.getB(coords).norm()
            else:
                if fieldCoord == Axis.X:
                    func = lambda self, coords: self.grid.getB(coords).x
                elif fieldCoord == Axis.Y:
                    func = lambda self, coords: self.grid.getB(coords).y
                elif fieldCoord == Axis.Z:
                    func = lambda self, coords: self.grid.getB(coords).z   
                    
        elif field == Field.J:
            if norm:
                func = lambda self, coords: self.grid.getJ(coords).norm()
            else:
                if fieldCoord == Axis.X:
                    func = lambda self, coords: self.grid.getJ(coords).x
                elif fieldCoord == Axis.Y:
                    func = lambda self, coords: self.grid.getJ(coords).y
                elif fieldCoord == Axis.Z:
                    func = lambda self, coords: self.grid.getJ(coords).z 
                    
        return func                    
    
    
    def getFieldPlane_(self, coords, shape, plane, lastCoordinateValue, func):
        field = np.zeros(shape=(shape[1], shape[0]))
        
        if plane == Plane.XOY:
            for i0 in range(shape[0]):
                for i1 in range(shape[1]):
                    coord = hichi.vector3d(coords[0][i0], coords[1][i1], lastCoordinateValue)
                    field[i1, i0] = func(self, coord) 
            return field
            
        elif plane == Plane.XOZ:
            for i0 in range(shape[0]):
                for i1 in range(shape[1]):
                    coord = hichi.vector3d(coords[0][i0], lastCoordinateValue, coords[1][i1])
                    field[i1, i0] = func(self, coord) 
            return field
            
        elif plane == Plane.YOZ:
            for i0 in range(shape[0]):
                for i1 in range(shape[1]):
                    coord = hichi.vector3d(lastCoordinateValue, coords[0][i0], coords[1][i1])
                    field[i1, i0] = func(self, coord) 
            return field
            
        return None
    
    
    def getFieldAxis_(self, coords, nPoints, axis, lastCoordinateValue, func):
        field = np.zeros(shape=(nPoints))
        
        if axis == Axis.X:
            for i in range(nPoints):
                coord = hichi.vector3d(coords[i], lastCoordinateValue[0], lastCoordinateValue[1])
                field[i] = func(self, coord) 
            return field
            
        elif axis == Axis.Y:
            for i in range(nPoints):
                coord = hichi.vector3d(lastCoordinateValue[0], coords[i], lastCoordinateValue[1])
                field[i] = func(self, coord) 
            return field
            
        elif axis == Axis.Z:
            for i in range(nPoints):
                coord = hichi.vector3d(lastCoordinateValue[0], lastCoordinateValue[1], coords[i])
                field[i] = func(self, coord) 
            return field
            
        return None


class Reader:
    
    def __init__(self, dir="./"):
        self.dir = os.path.join(os.getcwd(), dir)
        
    
    def readFile2d(self, nameFile="field.csv"):         
        with open(os.path.join(self.dir, nameFile), "r") as file:
            field = None
            lines = file.readlines()
            ny = len(lines)
            if ny > 1:
                nx = len(lines[0].split(";"))
                if nx > 1:
                    field = np.zeros(shape=(nx, ny))
                    for j in range(ny-1):
                        arr = lines[j].split(";") 
                        for i in range(nx-1):
                            field[i, j] = float(arr[i])
        
        return field
                  
        
    def readFile1d(self, nameFile="field.csv"):
        with open(os.path.join(self.dir, nameFile), "r") as file:
            field = None
            lines = file.readlines()
            nx = len(lines)
            if nx > 1:
                field = np.zeros(shape=(nx))
                for i in range(nx-1):
                    field[i] = float(lines[i].split()[0])
        
        return field     