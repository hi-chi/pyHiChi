import pyHiChi as hichi
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import os
from hichi_primitives import Axis, Plane, Field, getCoordValue   
    
    
class Visual:
    
    def __init__(self, grid, minCoords, maxCoords, dir="./", dpi=500, fontsize=17):
        self.grid = grid
        self.minCoords = minCoords
        self.maxCoords = maxCoords
        self.dir = os.path.join(os.getcwd(), dir)
        self.dpi = 500
        matplotlib.rcParams.update({"font.size" : fontsize})
        
    
    def savePictureInPlane(self, shape, plane=Plane.XOY, lastCoordinateValue=0.0,
                           field=Field.E, fieldCoord=Axis.X, norm=False,
                           valueLimits=(None, None),
                           namePicture="field.png"):
        title = ("$|%s|$" % (field.value)) if norm else ("$%s%s$" %(field.value, fieldCoord.value))
        xlabel = plane.value[0].value
        ylabel = plane.value[1].value
        minCoords = (getCoordValue(self.minCoords, plane.value[0]), getCoordValue(self.minCoords, plane.value[1]))
        maxCoords = (getCoordValue(self.maxCoords, plane.value[0]), getCoordValue(self.maxCoords, plane.value[1]))
        
        fig = plt.figure()
        ax, im = self.createAxPlane_(fig, shape, title, xlabel, ylabel, minCoords, maxCoords, valueLimits)
        
        coords0 = np.linspace(minCoords[0], maxCoords[0], shape[0])
        coords1 = np.linspace(minCoords[1], maxCoords[1], shape[1])         
        fields = self.getFieldPlane_((coords0, coords1), shape, plane, lastCoordinateValue,
            self.generateFuncGet_(field, fieldCoord, norm))
            
        im.set_array(fields)
        
        fig.tight_layout()
        
        plt.savefig(os.path.join(self.dir, namePicture), dpi=self.dpi)
        plt.close(fig=fig)
        
        
    def savePictureInAxis(self, nPoints, axis=Axis.X, lastCoordinateValue=(0.0, 0.0),
                          field=Field.E, fieldCoord=Axis.X, norm=False, linePlot="-", label="",
                          yLimits=None,
                          namePicture="field.png"):
        ylabel = ("$|%s|$" % (field.value)) if norm else ("$%s%s$" %(field.value, fieldCoord.value))
        xlabel = axis.value
        title = ""
        minCoords = getCoordValue(self.minCoords, axis)
        maxCoords = getCoordValue(self.maxCoords, axis)
        
        fig = plt.figure()
        ax = self.createAxAxis_(fig, title, xlabel, ylabel, minCoords, maxCoords, yLimits)
        
        coords = np.linspace(minCoords, maxCoords, nPoints)
        fields = self.getFieldAxis_(coords, nPoints, axis, lastCoordinateValue,
            self.generateFuncGet_(field, fieldCoord, norm))
            
        ax.plot(coords, fields, linePlot, label=label)
        
        fig.tight_layout()
                
        plt.savefig(os.path.join(self.dir, namePicture), dpi=self.dpi)
        plt.close(fig=fig)     


    def animateInPlane(self, funcUpdate, nIter, shape, plane=Plane.XOY, lastCoordinateValue=0.0,
                       field=Field.E, fieldCoord=Axis.X, norm=False,
                       valueLimits=(None, None), interval=1):
        title = ("$|%s|$" % (field.value)) if norm else ("$%s%s$" %(field.value, fieldCoord.value))
        xlabel = plane.value[0].value
        ylabel = plane.value[1].value
        minCoords = (getCoordValue(self.minCoords, plane.value[0]), getCoordValue(self.minCoords, plane.value[1]))
        maxCoords = (getCoordValue(self.maxCoords, plane.value[0]), getCoordValue(self.maxCoords, plane.value[1]))
        
        fig = plt.figure()
        ax, im = self.createAxPlane_(fig, shape, title, xlabel, ylabel, minCoords, maxCoords, valueLimits)
        
        coords0 = np.linspace(minCoords[0], maxCoords[0], shape[0])
        coords1 = np.linspace(minCoords[1], maxCoords[1], shape[1])         

        fig.tight_layout()
                
        def animate_(i):
            if (i > nIter):
                exit()
            funcUpdate()
            fields = self.getFieldPlane_((coords0, coords1), shape, plane, lastCoordinateValue,
                self.generateFuncGet_(field, fieldCoord, norm))
            im.set_array(fields)
            return im,   
    
        ani = animation.FuncAnimation(fig, animate_, interval=interval, blit=True)
        plt.show()  


    def animateInAxis(self, funcUpdate, nIter, nPoints, axis=Axis.X, lastCoordinateValue=(0.0, 0.0),
                      field=Field.E, fieldCoord=Axis.X, norm=False, linePlot="-", label="",
                      yLimits=None, interval=10):
        ylabel = ("$|%s|$" % (field.value)) if norm else ("$%s%s$" %(field.value, fieldCoord.value))
        xlabel = axis.value
        title = ""
        minCoords = getCoordValue(self.minCoords, axis)
        maxCoords = getCoordValue(self.maxCoords, axis)
        
        fig = plt.figure()
        ax = self.createAxAxis_(fig, title, xlabel, ylabel, minCoords, maxCoords, yLimits)
        
        coords = np.linspace(minCoords, maxCoords, nPoints)
        fields = self.getFieldAxis_(coords, nPoints, axis, lastCoordinateValue,
            self.generateFuncGet_(field, fieldCoord, norm))
            
        line, = ax.plot(coords, fields, linePlot, label=label)
        
        fig.tight_layout()               
                
        def animate_(i):
            if (i > nIter):
                exit()
            funcUpdate()
            fields = self.getFieldAxis_(coords, nPoints, axis, lastCoordinateValue,
                self.generateFuncGet_(field, fieldCoord, norm))
            line.set_data(coords, fields)
            return line,   
    
        ani = animation.FuncAnimation(fig, animate_, interval=interval, blit=True)
        plt.show()   
        
        
    def createAxPlane_(self, fig, shape, title, xlabel, ylabel, minCoords, maxCoords, valueLimits):
        fig.clear()
        ax = fig.add_subplot(1, 1, 1)
        ax.title.set_text(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.tick_params(axis='both', which='major')
        im = ax.imshow(np.zeros(shape=shape), cmap='RdBu', interpolation='none',
            extent=(minCoords[0], maxCoords[0], minCoords[1], maxCoords[1]),\
            animated = True, aspect='auto', vmin=valueLimits[0], vmax=valueLimits[1])
        fig.colorbar(im, ax=ax)
        return ax, im    
        
        
    def createAxAxis_(self, fig, title, xlabel, ylabel, minCoords, maxCoords, yLimits):
        fig.clear()
        ax = fig.add_subplot(1, 1, 1)
        ax.title.set_text(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.grid()
        ax.set_xlim(minCoords, maxCoords)
        if yLimits: ax.set_ylim(yLimits)
        ax.tick_params(axis='both', which='major')
        return ax        
        
        
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
                    field[shape[1]-i1-1, i0] = func(self, coord) 
            return field
            
        elif plane == Plane.XOZ:
            for i0 in range(shape[0]):
                for i1 in range(shape[1]):
                    coord = hichi.vector3d(coords[0][i0], lastCoordinateValue, coords[1][i1])
                    field[shape[1]-i1-1, i0] = func(self, coord) 
            return field
            
        elif plane == Plane.YOZ:
            for i0 in range(shape[0]):
                for i1 in range(shape[1]):
                    coord = hichi.vector3d(lastCoordinateValue, coords[0][i0], coords[1][i1])
                    field[shape[1]-i1-1, i0] = func(self, coord) 
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
    
    