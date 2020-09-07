import pyHiChi as hichi
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import os
from hichi_primitives import Axis, Plane, Field, get_coord_value   
    
    
class Visual:
    
    def __init__(self, field, min_coords, max_coords, dir="./", dpi=500, fontsize=17):
        self.field = field
        self.min_coords = min_coords
        self.max_coords = max_coords
        self.dir = os.path.join(os.getcwd(), dir)
        self.dpi = 500
        matplotlib.rcParams.update({"font.size" : fontsize})
        
    
    def save_plane_to_image(self, shape, plane=Plane.XOY, last_coordinate_value=0.0,
                            field=Field.E, field_coord=Axis.X, norm=False,
                            value_limits=(None, None),
                            name_picture="field.png"
                           ):
        title = ("$|%s|$" % (field.value)) if norm else ("$%s%s$" %(field.value, field_coord.value))
        xlabel = plane.value[0].value
        ylabel = plane.value[1].value
        min_coords = (get_coord_value(self.min_coords, plane.value[0]), get_coord_value(self.min_coords, plane.value[1]))
        max_coords = (get_coord_value(self.max_coords, plane.value[0]), get_coord_value(self.max_coords, plane.value[1]))
        
        fig = plt.figure()
        ax, im = self.create_ax_plane_(fig, shape, title, xlabel, ylabel, min_coords, max_coords, value_limits)
        
        coords0 = np.linspace(min_coords[0], max_coords[0], shape[0])
        coords1 = np.linspace(min_coords[1], max_coords[1], shape[1])         
        fields = self.get_field_plane_((coords0, coords1), shape, plane, last_coordinate_value,
            self.generate_func_get_(field, field_coord, norm))
            
        im.set_array(fields)
        
        fig.tight_layout()
        
        plt.savefig(os.path.join(self.dir, name_picture), dpi=self.dpi)
        plt.close(fig=fig)
        
        
    def save_axis_to_image(self, n_points, axis=Axis.X, last_coordinate_value=(0.0, 0.0),
                           field=Field.E, field_coord=Axis.X, norm=False, line_plot="-", label="",
                           y_limits=None,
                           name_picture="field.png"
                          ):
        ylabel = ("$|%s|$" % (field.value)) if norm else ("$%s%s$" %(field.value, field_coord.value))
        xlabel = axis.value
        title = ""
        min_coords = get_coord_value(self.min_coords, axis)
        max_coords = get_coord_value(self.max_coords, axis)
        
        fig = plt.figure()
        ax = self.create_ax_axis_(fig, title, xlabel, ylabel, min_coords, max_coords, y_limits)
        
        coords = np.linspace(min_coords, max_coords, n_points)
        fields = self.get_field_axis_(coords, n_points, axis, last_coordinate_value,
            self.generate_func_get_(field, field_coord, norm))
            
        ax.plot(coords, fields, line_plot, label=label)
        
        fig.tight_layout()
                
        plt.savefig(os.path.join(self.dir, name_picture), dpi=self.dpi)
        plt.close(fig=fig)     


    def animate_plane(self, func_update, n_iter, shape, plane=Plane.XOY, last_coordinate_value=0.0,
                      field=Field.E, field_coord=Axis.X, norm=False,
                      value_limits=(None, None), interval=1
                     ):
        title = ("$|%s|$" % (field.value)) if norm else ("$%s%s$" %(field.value, field_coord.value))
        xlabel = plane.value[0].value
        ylabel = plane.value[1].value
        min_coords = (get_coord_value(self.min_coords, plane.value[0]), get_coord_value(self.min_coords, plane.value[1]))
        max_coords = (get_coord_value(self.max_coords, plane.value[0]), get_coord_value(self.max_coords, plane.value[1]))
        
        fig = plt.figure()
        ax, im = self.create_ax_plane_(fig, shape, title, xlabel, ylabel, min_coords, max_coords, value_limits)
        
        coords0 = np.linspace(min_coords[0], max_coords[0], shape[0])
        coords1 = np.linspace(min_coords[1], max_coords[1], shape[1])         

        fig.tight_layout()
                
        def animate_(i):
            if (i > n_iter):
                exit()
            func_update()
            fields = self.get_field_plane_((coords0, coords1), shape, plane, last_coordinate_value,
                self.generate_func_get_(field, field_coord, norm))
            im.set_array(fields)
            return im,   
    
        ani = animation.FuncAnimation(fig, animate_, interval=interval, blit=True)
        plt.show()  


    def animate_axis(self, func_update, n_iter, n_points, axis=Axis.X, last_coordinate_value=(0.0, 0.0),
                     field=Field.E, field_coord=Axis.X, norm=False, line_plot="-", label="",
                     y_limits=None, interval=10
                    ):
        ylabel = ("$|%s|$" % (field.value)) if norm else ("$%s%s$" %(field.value, field_coord.value))
        xlabel = axis.value
        title = ""
        min_coords = get_coord_value(self.min_coords, axis)
        max_coords = get_coord_value(self.max_coords, axis)
        
        fig = plt.figure()
        ax = self.create_ax_axis_(fig, title, xlabel, ylabel, min_coords, max_coords, y_limits)
        
        coords = np.linspace(min_coords, max_coords, n_points)
        fields = self.get_field_axis_(coords, n_points, axis, last_coordinate_value,
            self.generate_func_get_(field, field_coord, norm))
            
        line, = ax.plot(coords, fields, line_plot, label=label)
        
        fig.tight_layout()               
                
        def animate_(i):
            if (i > n_iter):
                exit()
            func_update()
            fields = self.get_field_axis_(coords, n_points, axis, last_coordinate_value,
                self.generate_func_get_(field, field_coord, norm))
            line.set_data(coords, fields)
            return line,   
    
        ani = animation.FuncAnimation(fig, animate_, interval=interval, blit=True)
        plt.show()   
        
        
    def create_ax_plane_(self, fig, shape, title, xlabel, ylabel, min_coords, max_coords, value_limits):
        fig.clear()
        ax = fig.add_subplot(1, 1, 1)
        ax.title.set_text(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.tick_params(axis='both', which='major')
        im = ax.imshow(np.zeros(shape=shape), cmap='RdBu', interpolation='none',
            extent=(min_coords[0], max_coords[0], min_coords[1], max_coords[1]),\
            animated = True, aspect='auto', vmin=value_limits[0], vmax=value_limits[1])
        fig.colorbar(im, ax=ax)
        return ax, im    
        
        
    def create_ax_axis_(self, fig, title, xlabel, ylabel, min_coords, max_coords, y_limits):
        fig.clear()
        ax = fig.add_subplot(1, 1, 1)
        ax.title.set_text(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.grid()
        ax.set_xlim(min_coords, max_coords)
        if y_limits: ax.set_ylim(y_limits)
        ax.tick_params(axis='both', which='major')
        return ax        
        
        
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
                    field[shape[1]-i1-1, i0] = func(self, coord) 
            return field
            
        elif plane == Plane.XOZ:
            for i0 in range(shape[0]):
                for i1 in range(shape[1]):
                    coord = hichi.Vector3d(coords[0][i0], last_coordinate_value, coords[1][i1])
                    field[shape[1]-i1-1, i0] = func(self, coord) 
            return field
            
        elif plane == Plane.YOZ:
            for i0 in range(shape[0]):
                for i1 in range(shape[1]):
                    coord = hichi.Vector3d(last_coordinate_value, coords[0][i0], coords[1][i1])
                    field[shape[1]-i1-1, i0] = func(self, coord) 
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
    
    