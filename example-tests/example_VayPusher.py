import sys
sys.path.append("../bin/")
import pyHiChi as hichi
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as animation

#Particle initialization

particle = hichi.Particle()
pos = hichi.Vector3d(0, 0, 0)
V = hichi.Vector3d(hichi.LIGHT_VELOCITY *10**-2, 0, 0)
particle.set_velocity(V)
particle.set_position(pos)
p0=particle.get_momentum()
x_coords, y_coords = [0], [0]

E = hichi.Vector3d(0, 0, 0)
B = hichi.Vector3d(0, 0, 0.1)
field = hichi.FieldValue(E, B)
N=100
timeStep = np.pi*particle.get_mass()*hichi.LIGHT_VELOCITY*particle.get_gamma()/(abs(particle.get_charge()* 0.1 * N))

# Сalculation by the Vay method

pusher = hichi.VayPusher()
for i in range(2*N):
    pusher(particle, field, timeStep)
    pos_real = particle.get_position()
    #plt.plot(pos_real.x,pos_real.y,'ro')
    x_coords.append(pos_real.x)
    y_coords.append(pos_real.y)

# Analytical coordinates of the particle

pos_y = -2 * p0.x * hichi.LIGHT_VELOCITY/(particle.get_charge()* 0.1)
pos_analytical=hichi.Vector3d(0, pos_y, 0)

# Animation

fig = plt.figure(figsize=(9, 6))
ax = plt.axes(xlim=(-250, 250), ylim=(-25, 375))
line, = ax.plot([], [], lw=2, color ='r')

t = np.arange(0, 2*np.pi, 0.01)
r = pos_y/2
plt.plot(r*np.sin(t), pos_y/2 + r*np.cos(t), lw=2 , color='b')
plt.axis('scaled')
 
def init():
    line.set_data([], [])
    return line,
 
xdata, ydata = [], []
 
def animate(i):
    t = 0.1 * i
    xdata.append(x_coords[i])
    ydata.append(y_coords[i])
    line.set_data(xdata, ydata)
    return line,
 
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=2*N, interval=30, blit=True)

plt.legend(['Real path', 'Analytical path'], loc='lower center', bbox_to_anchor=(1.2, 0.5))

plt.show()