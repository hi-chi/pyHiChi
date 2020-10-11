import sys
sys.path.append("../../bin")
sys.path.append("../../python_modules")

import plotly
import plotly.figure_factory as ff
import plotly.graph_objects as go
import plotly.express as px
import numpy as np
from scipy.spatial import Delaunay

import pyHiChi as hichi
from tight_focusing_fields import SphericalPulsePython

import math

spherical_pulse = SphericalPulsePython(f_number = 0.25,
                                 R0 = 16,
                                 pulselength = 2.0,
                                 edge_smoothing_angle = 0.0
                                )
MAX_THETA = spherical_pulse.opening_angle
R0 = spherical_pulse.R0


def get_r_by_y(y, R):
    return (R**2-y**2)**0.5

def get_theta_by_y(y, R):
    v = np.arccos(abs((R*np.cos(MAX_THETA))/(get_r_by_y(y, R))))
    return v if not math.isnan(v) else 0.0

def get_phi_of_0_by_R(R):
    return np.arcsin(np.abs(0/(R*np.sin(MAX_THETA))))

R_MAX = spherical_pulse.R0 + 0.5*spherical_pulse.pulselength
R_MIN = spherical_pulse.R0 - 0.5*spherical_pulse.pulselength
D = 2.3*spherical_pulse.pulselength

mapping = hichi.TightFocusingMapping(spherical_pulse.R0, spherical_pulse.pulselength, D)
X_MIN = mapping.get_min_coord()
X_MAX = mapping.get_max_coord()

MIN_COORDS = hichi.Vector3d(-25, -20, -20)
MAX_COORDS = hichi.Vector3d(1, 20, 20)

FACTOR = 1 # 10

main_colormap = px.colors.sequential.YlOrBr_r
secondary_colormap = px.colors.sequential.PuBu_r
gray_colormap = px.colors.sequential.Greys_r[1:len(px.colors.sequential.Greys_r)-2]
purple_colormap = px.colors.sequential.Purples_r[1:len(px.colors.sequential.Purples_r)-2]


fig_arr = []
aspectratio=dict(x=1,
                 y=(MAX_COORDS.y-MIN_COORDS.y)/(MAX_COORDS.x-MIN_COORDS.x),
                 z=(MAX_COORDS.y-MIN_COORDS.y)/(MAX_COORDS.x-MIN_COORDS.x)
                )
                                
# ------------------- define area ---------------------


def define_area():
    x = [MIN_COORDS.x, MIN_COORDS.x, MIN_COORDS.x, MIN_COORDS.x, MAX_COORDS.x, MAX_COORDS.x, MAX_COORDS.x, MAX_COORDS.x]
    y = [MIN_COORDS.y, MIN_COORDS.y, MAX_COORDS.y, MAX_COORDS.y, MIN_COORDS.y, MIN_COORDS.y, MAX_COORDS.y, MAX_COORDS.y]
    z = [MIN_COORDS.z, MAX_COORDS.z, MIN_COORDS.z, MAX_COORDS.z, MIN_COORDS.z, MAX_COORDS.z, MIN_COORDS.z, MAX_COORDS.z]
    fig = go.Figure(go.Scatter3d(x=x, y=y, z=z,
                                 showlegend=False,
                                )
                   )
    fig.data[0].update(opacity=0.0)
    fig_arr.append(fig)


# ------------------- compute -------------------------

def create_2d_points(u_min, u_max, Nu, v_min_f, v_max_f, Nv):
    u_arr = np.linspace(u_min, u_max, Nu)
    u = np.array([[u for i in range(Nv)] for u in u_arr])
    v = np.array([np.linspace(v_min_f(u), v_max_f(u), Nv) for u in u_arr])
    u = u.flatten()
    v = v.flatten()
    return u, v


def get_points_of_sphere(R):
    y, theta = create_2d_points(-R*np.sin(MAX_THETA), R*np.sin(MAX_THETA), 30*FACTOR,
                                lambda y: 0 if y<0 else -get_theta_by_y(y, R),
                                lambda y: get_theta_by_y(y, R),
                                30*FACTOR
                               )
    
    x = -get_r_by_y(y, R)*np.cos(theta)
    z = -get_r_by_y(y, R)*np.sin(theta)
    
    return y, theta, x, y, z


def get_points_of_edge():
    min_v_f = lambda R: 0
    max_v_f = lambda R: 3/2*np.pi
    R, phi = create_2d_points(R_MIN, R_MAX, 3*FACTOR,
                              min_v_f, max_v_f, 30*FACTOR
                             )
    phi = phi[::-1]
    
    x = -R*np.cos(MAX_THETA)
    y = -R*np.cos(phi)*np.sin(MAX_THETA)
    z = -R*np.sin(phi)*np.sin(MAX_THETA)
    
    return R, phi, x, y, z
    
    
def get_points_of_section(axis):
    R, r = create_2d_points(R_MIN, R_MAX, 30*FACTOR,
                            lambda R: get_r_by_y(0, R),
                            lambda R: get_r_by_y(0, R),
                            30*FACTOR
                           )

    R, theta = create_2d_points(R_MIN, R_MAX, 30*FACTOR,
                                lambda R: -get_theta_by_y(0, R),
                                lambda R: 0,
                                30*FACTOR
                               )
                               
    if axis=="y":                           
        x = -r*np.cos(theta)
        y = np.full(x.shape, 0)
        z = -r*np.sin(theta)
    else:
        x = -r*np.cos(theta)
        z = np.full(x.shape, 0)
        y = r*np.sin(theta)
    
    return r, theta, x, y, z


def get_triangles(u, v):
    points2D = np.vstack([u, v]).T
    tri = Delaunay(points2D)
    simplices = tri.simplices
    return tri.simplices


def combine_points(data):
    x = np.concatenate([s[2] for s in data])
    y = np.concatenate([s[3] for s in data])
    z = np.concatenate([s[4] for s in data])
    
    simplices = np.ndarray(shape=(0,3)).astype(int)
    for i in range(len(data)):
        tmp_tri = get_triangles(data[i][0], data[i][1])
        j = 0
        while j < len(tmp_tri):
            p0 = [data[i][2][tmp_tri[j][0]], data[i][3][tmp_tri[j][0]], data[i][4][tmp_tri[j][0]]]
            p1 = [data[i][2][tmp_tri[j][1]], data[i][3][tmp_tri[j][1]], data[i][4][tmp_tri[j][1]]]
            p2 = [data[i][2][tmp_tri[j][2]], data[i][3][tmp_tri[j][2]], data[i][4][tmp_tri[j][2]]]
            if ((p0[0]-p1[0])**2+(p0[1]-p1[1])**2+(p0[2]-p1[2])**2 > 2**2) or\
                ((p1[0]-p2[0])**2+(p1[1]-p2[1])**2+(p1[2]-p2[2])**2 > 2**2) or\
                ((p0[0]-p2[0])**2+(p0[1]-p2[1])**2+(p0[2]-p2[2])**2 > 2**2):
                tmp_tri = np.delete(tmp_tri, j, axis=0)
            else: j += 1
        tmp_tri = tmp_tri + np.full(tmp_tri.shape, sum((len(data[j][2]) for j in range(i))))
        simplices = np.concatenate([simplices, tmp_tri])
    
    return x, y, z, simplices
    

def create_pulse(shift):  
    s1 = get_points_of_sphere(R_MIN)
    s2 = get_points_of_edge()
    s3 = get_points_of_sphere(R_MAX)
    
    x, y, z, simplices = combine_points([s1, s2, s3])
    x = x + np.full(x.shape, shift)
    
    return x, y, z, simplices
    

def create_section(shift):
    s1 = get_points_of_section("y")
    s2 = get_points_of_section("z")
    
    x, y, z, simplices = combine_points([s1, s2])
    x = x + np.full(x.shape, shift)
    
    return x, y, z, simplices
    

# -------------------------- plot ---------------------------------


def draw_circle(x, r, line):
    N = 30*FACTOR
    phi = np.linspace(0, np.pi/2, N)
    x = np.full((N), x)
    y = -r*np.sin(phi)
    z = -r*np.cos(phi)
    fig = go.Figure(go.Scatter3d(x=x, y=y, z=z,
                                 mode="lines",
                                 line=line,
                                 showlegend=False,
                                )
                   )
    fig_arr.append(fig)


def draw_plane(axis, points, colormap, opacity, if_show_circle=False):
    axis1, axis2 = (0, 1) if axis==2 else (0, 2) if axis == 1 else (1, 2)
    N = 50*FACTOR
    u = np.linspace(min([p[axis1] for p in points]), max([p[axis1] for p in points]), N)
    v = np.linspace(min([p[axis2] for p in points]), max([p[axis2] for p in points]), N)
    u, v = np.meshgrid(u, v)
    u = list(u.flatten())
    v = list(v.flatten())
    
    c = points[0][axis]
    x = list(u) if axis1==0 else list(v) if axis2==0 else list(np.full(len(u), c))
    y = list(u) if axis1==1 else list(v) if axis2==1 else list(np.full(len(u), c))
    z = list(u) if axis1==2 else list(v) if axis2==2 else list(np.full(len(u), c))
    
    i = 0
    while i <len(x):
        if y[i]<0 and z[i]>0:
            u.pop(i)
            v.pop(i)
            x.pop(i)
            y.pop(i)
            z.pop(i)
        else: i += 1
    
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    u = np.array(u)
    v = np.array(v)
    
    simplices = get_triangles(u, v)
    
    fig = ff.create_trisurf(x=x, y=y, z=z,
                            simplices=simplices,
                            aspectratio=aspectratio,
                            plot_edges=False,
                            colormap=colormap[3:],
                            show_colorbar=False,
                           )
    fig.data[0].update(opacity=opacity)
    fig_arr.append(fig)
    
    line = dict(width=2, color="black")
    fig = go.Figure(go.Scatter3d(x=[p[0] for p in points]+[points[0][0]],
                                 y=[p[1] for p in points]+[points[0][1]],
                                 z=[p[2] for p in points]+[points[0][2]],
                                 mode="lines",
                                 line=line,
                                 showlegend=False,
                                )
                   )
    fig_arr.append(fig)
    
    if (if_show_circle):
        line = dict(width=2, color="black", dash="longdash")
        R = R_MAX + 0.02
        
        def plot_circle(shift):
            theta = np.arccos(abs((X_MAX + shift)/R))
            x = -R*np.cos(theta)-(D if c == X_MIN else 0)
            r = -R*np.sin(theta)
            draw_circle(x-shift, r, line)
        
        plot_circle(0)
        plot_circle(D)
        

def draw_edge_line(R, shift):
    N = 30*FACTOR
    phi = np.linspace(0, 3/2*np.pi, N)
    x = np.full((N), -R*np.cos(MAX_THETA)+shift)
    y = -R*np.sin(MAX_THETA)*np.cos(phi)
    z = -R*np.sin(MAX_THETA)*np.sin(phi)

    fig = go.Figure(go.Scatter3d(x=x, y=y, z=z,
                                 mode="lines",
                                 line=dict(width=2, color="rgb(100, 100, 100)"),
                                 showlegend=False,
                                )
                   )
    if shift == D: fig.data[0].update(opacity=0.4)
    fig_arr.append(fig)

    
def draw_sphere(shift, colormap):
    x, y, z, simplices = create_pulse(shift)
    
    fig = ff.create_trisurf(x=x, y=y, z=z,
                            simplices=simplices,
                            aspectratio=aspectratio,
                            plot_edges=False,
                            colormap=colormap,
                            show_colorbar=False,
                           )
    if shift == D: fig.data[0].update(opacity=0.4)
    fig_arr.append(fig)


def draw_section(shift, colormap):
    x, y, z, simplices = create_section(shift)
    
    fig = ff.create_trisurf(x=x, y=y, z=z,
                            simplices=simplices,
                            aspectratio=aspectratio,
                            plot_edges=False,
                            colormap=colormap,
                            show_colorbar=False,
                            color_func=lambda x, y, z: np.abs(spherical_pulse.mask(x-shift, y, z)),
                           )
    if shift == D: fig.data[0].update(opacity=0.4)
    fig_arr.append(fig)


def draw_pulse(shift, colormap):
    draw_sphere(shift, colormap)
    draw_section(shift, colormap)
    draw_edge_line(R_MIN, shift)
    draw_edge_line(R_MAX, shift)
    
    
annotations = []
def draw_label(x1, x2, y, z, text, xlabel=None, ylabel=None):
    fig = go.Figure(go.Scatter3d(
        x=[x1, x2],
        y=[y, y],
        z=[z, z],
        mode="lines",
        line=dict(color="rgb(100,100,100)"),
        showlegend=False,
    ))
    if not xlabel: xlabel = (x2+x1)/2
    if not ylabel: ylabel = y
    annotations.append(
            dict(x=xlabel, y=y-5, z=z+2,
                 showarrow=False,
                 text=text,
                 font=dict(
                     color="black",
                     size=12
                 ),
            )
        )
    fig_arr.append(fig)

def draw_label_pulses(x, y, z, text, ax=0, ay=-50):
    annotations.append(dict(
         x=x,
         y=y,
         z=z,
         ax=ax,
         ay=ay,
         font=dict(color="black",
                   size=12
                  ),
         arrowwidth=1,
         text=text,
         #arrowhead=1,
         textangle=0,
        ))
        

def add_focus():
    fig = go.Figure(go.Scatter3d(
        x=[0],
        y=[0],
        z=[0],
        mode="markers",
        showlegend=False,
        marker=dict(color="black", size=3),
    ))
    fig_arr.append(fig)

# ------------------------- run -----------------------------


define_area()

shift_arr = [i*D for i in range(-2, 1)]
colormap_arr = [main_colormap if s==0 else secondary_colormap for s in shift_arr]

for shift, colormap in zip(shift_arr, colormap_arr):
    draw_pulse(shift, colormap)

def get_plane_x(x, c):   
    return [[x, -c, -c],
              [x, -c, 0],
              [x, 0, c],
              [x, c, c],
              [x, c, -c],
             ]
           
draw_plane(0, get_plane_x(X_MIN, 18), purple_colormap, 0.6, False)
draw_plane(0, get_plane_x(X_MAX, 18), purple_colormap, 0.6, False)
# draw_plane(0, get_plane_x(-20, 17), gray_colormap, 0.5)
# draw_plane(0, get_plane_x(20, 17), gray_colormap, 0.5)

# draw_label(-20, 20, -13, -17, "full area")
#draw_label(X_MIN, X_MAX, -18, -18, "periodic subregion", -17, -20)

#draw_label_pulses(-8, 2, 16, "initial pulse", 10, -50)
#draw_label_pulses(-12, 6, 15, "repeated pulses", 0, -70)
#draw_label_pulses(-16, 6, 15, "", 60, -65)

add_focus()

fig = go.Figure()
fig.add_traces([f.data[0] for f in fig_arr])

for f in fig_arr:
    fig.layout.update(f.layout)
    
fig.update_layout(
    scene=dict(
        yaxis = dict(
            tickmode = 'array',
            tickvals = [15, 10, 5, 0, -5, -10, -15, -20],
        ),
        annotations=annotations,
        camera=dict(
            eye=dict(
                x=-0.6,
                y=-3.6,
                z=1.5,
            ),
        ),
    ),
    title={'text': ""}
)

fig.show()
#fig.write_image("3d_image.svg")

