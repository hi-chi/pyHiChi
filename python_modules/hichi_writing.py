import pyHiChi as hichi
import numpy as np
import os
from hichi_primitives import Axis, Plane, Field, get_coord_value


class Writer:
    
    def __init__(self, field, min_coords, max_coords, dir="./"):
        self.field = field
        self.min_coords = min_coords
        self.max_coords = max_coords
        self.dir = os.path.join(os.getcwd(), dir)
           
    def save_plane_to_file(self, shape, plane=Plane.XOY, last_coordinate_value=0.0,
                           field=Field.E, field_coord=Axis.X, norm=False,
                           name_file="field.csv"
                          ):
        min_coords = (get_coord_value(self.min_coords, plane.value[0]), get_coord_value(self.min_coords, plane.value[1]))
        max_coords = (get_coord_value(self.max_coords, plane.value[0]), get_coord_value(self.max_coords, plane.value[1]))
        coords0 = np.linspace(min_coords[0], max_coords[0], shape[0])
        coords1 = np.linspace(min_coords[1], max_coords[1], shape[1])         
        fields = self.get_field_plane_((coords0, coords1), shape, plane, last_coordinate_value,
            self.generate_func_get_(field, field_coord, norm))
            
        with open(os.path.join(self.dir, name_file), "w") as file:
            for iy in range(shape[1]):
                for ix in range(shape[0]):
                    file.write("%f;" % fields[iy, ix])
                file.write("\n")
                   
    def save_axis_to_file(self, n_points, axis=Axis.X, last_coordinate_value=(0.0, 0.0),
                          field=Field.E, field_coord=Axis.X, norm=False,
                          name_file="field.csv"
                         ):
        min_coords = get_coord_value(self.min_coords, axis)
        max_coords = get_coord_value(self.max_coords, axis)
        coords = np.linspace(min_coords, max_coords, n_points)
        fields = self.get_field_axis_(coords, n_points, axis, last_coordinate_value,
            self.generate_func_get_(field, field_coord, norm))
            
        with open(os.path.join(self.dir, name_file), "w") as file:
            for ix in range(n_points):
                file.write("%f\n" % fields[ix])
                   
    def generate_func_get_(self, field, field_coord, norm):
        func = None
        
        if field == Field.E:
            if norm:
                func = lambda self, coords: self.field.get_E(coords).norm()
            else:
                if field_coord == Axis.X:
                    func = lambda self, coords: self.field.get_E(coords).x
                elif field_coord == Axis.Y:
                    func = lambda self, coords: self.field.get_E(coords).y
                elif field_coord == Axis.Z:
                    func = lambda self, coords: self.field.get_E(coords).z
                    
        elif field == Field.B:
            if norm:
                func = lambda self, coords: self.field.get_B(coords).norm()
            else:
                if field_coord == Axis.X:
                    func = lambda self, coords: self.field.get_B(coords).x
                elif field_coord == Axis.Y:
                    func = lambda self, coords: self.field.get_B(coords).y
                elif field_coord == Axis.Z:
                    func = lambda self, coords: self.field.get_B(coords).z   
                    
        elif field == Field.J:
            if norm:
                func = lambda self, coords: self.field.get_J(coords).norm()
            else:
                if field_coord == Axis.X:
                    func = lambda self, coords: self.field.get_J(coords).x
                elif field_coord == Axis.Y:
                    func = lambda self, coords: self.field.get_J(coords).y
                elif field_coord == Axis.Z:
                    func = lambda self, coords: self.field.get_J(coords).z 
                    
        return func                    
       
    def get_field_plane_(self, coords, shape, plane, last_coordinate_value, func):
        field = np.zeros(shape=(shape[1], shape[0]))
        
        if plane == Plane.XOY:
            for i0 in range(shape[0]):
                for i1 in range(shape[1]):
                    coord = hichi.Vector3d(coords[0][i0], coords[1][i1], last_coordinate_value)
                    field[i1, i0] = func(self, coord) 
            return field
            
        elif plane == Plane.XOZ:
            for i0 in range(shape[0]):
                for i1 in range(shape[1]):
                    coord = hichi.Vector3d(coords[0][i0], last_coordinate_value, coords[1][i1])
                    field[i1, i0] = func(self, coord) 
            return field
            
        elif plane == Plane.YOZ:
            for i0 in range(shape[0]):
                for i1 in range(shape[1]):
                    coord = hichi.Vector3d(last_coordinate_value, coords[0][i0], coords[1][i1])
                    field[i1, i0] = func(self, coord) 
            return field
            
        return None
       
    def get_field_axis_(self, coords, n_points, axis, last_coordinate_value, func):
        field = np.zeros(shape=(n_points))
        
        if axis == Axis.X:
            for i in range(n_points):
                coord = hichi.Vector3d(coords[i], last_coordinate_value[0], last_coordinate_value[1])
                field[i] = func(self, coord) 
            return field
            
        elif axis == Axis.Y:
            for i in range(n_points):
                coord = hichi.Vector3d(last_coordinate_value[0], coords[i], last_coordinate_value[1])
                field[i] = func(self, coord) 
            return field
            
        elif axis == Axis.Z:
            for i in range(n_points):
                coord = hichi.Vector3d(last_coordinate_value[0], last_coordinate_value[1], coords[i])
                field[i] = func(self, coord) 
            return field
            
        return None


class Reader:
    
    def __init__(self, dir="./"):
        self.dir = os.path.join(os.getcwd(), dir)
           
    def read_file_2d(self, name_file="field.csv"):         
        with open(os.path.join(self.dir, name_file), "r") as file:
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
                        
    def read_file_1d(self, name_file="field.csv"):
        with open(os.path.join(self.dir, name_file), "r") as file:
            field = None
            lines = file.readlines()
            nx = len(lines)
            if nx > 1:
                field = np.zeros(shape=(nx))
                for i in range(nx-1):
                    field[i] = float(lines[i].split()[0])
        
        return field
        