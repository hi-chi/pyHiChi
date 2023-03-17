import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys
import re


LIGHT_VELOCITY = 29979245800.0
TIME_STEP = 0.1 / LIGHT_VELOCITY  # time step to draw animation, s

# first arg is file with geant output
if len(sys.argv) < 2:
    print("Expected file name as command line argument")
    exit()

# geant track table 
class ParticleInfo:
    def __init__(self):
        self.picID = 0
        self.parentID = 0
        self.pic_type = ""
        self.position_x = 0.               # mm
        self.position_y = 0.               # mm
        self.position_z = 0.               # mm
        self.global_time = 0.              # ns
        self.momentum_direct_x = 0.
        self.momentum_direct_y = 0.
        self.momentum_direct_z = 0.
        self.momentum_x = 0.               # Mev/c
        self.momentum_y = 0.               # Mev/c
        self.momentum_z = 0.               # Mev/c
        self.velocity = 0.                 # mm/ns
        self.weight = 0.
    
    def convert_units(self):
        self.position_x *= 1e-1            # cm
        self.position_y *= 1e-1            # cm
        self.position_z *= 1e-1            # cm
        self.global_time *= 1e-9           # s
        self.momentum_x *= 1.60218e-6      # erg/c?
        self.momentum_y *= 1.60218e-6      # erg/c?
        self.momentum_z *= 1.60218e-6      # erg/c?
        self.velocity *= 1e-1/1e-9         # cm/s
        
    def __str__(self):
        return "(t=%e, x=%0.3f, y=%0.3f, z=%0.3f)" % \
            (self.global_time, self.position_x, self.position_y, self.position_z)
    def __repr__(self):
        return str(self)
              

# track point to draw
class TrackPoint:
    def __init__(self, time, x, y, z):   
        self.time = time
        self.x = x
        self.y = y
        self.z = z
        
    def __str__(self):
        return "(t=%e, x=%0.3f, y=%0.3f, z=%0.3f)" % \
            (self.time, self.x, self.y, self.z)
    def __repr__(self):
        return str(self)     


# keeps full geant track as array of ParticleTrackStep
# keeps full track to draw as array of TrackPoint
class ParticleTrack:
    def __init__(self):
        self.geant_steps = []
        self.track = []
    
    def shift_time(self, time):
        for step in self.geant_steps:
            step.global_time -= time
    
    # converts geant track to draw track
    # draw track keeps more points than geant track
    def compute_track(self, TIME_STEP, t_start, t_final):
        start = self.geant_steps[0]
        
        t_cur = t_start
        for t in np.arange(t_start, start.global_time, TIME_STEP):
            self.track.append(TrackPoint(t, float("inf"), float("inf"), float("inf")))
            t_cur = t
        t_cur += TIME_STEP
        x = start.position_x + start.velocity * start.momentum_direct_x * (t_cur - start.global_time)
        y = start.position_y + start.velocity * start.momentum_direct_y * (t_cur - start.global_time)
        z = start.position_z + start.velocity * start.momentum_direct_z * (t_cur - start.global_time)
        self.track.append(TrackPoint(t_cur, x, y, z))
        
        index_cur_step = 0
        current_step = self.geant_steps[index_cur_step].final_pi
        
        for t in np.arange(self.track[-1].time, t_final, TIME_STEP):
            if index_cur_step + 1 < len(self.geant_steps) and t > self.geant_steps[index_cur_step].final_pi.global_time:
                index_cur_step += 1
                current_step = self.geant_steps[index_cur_step]
                x = current_step.position_x + current_step.velocity * current_step.momentum_direct_x * (t - current_step.global_time)
                y = current_step.position_y + current_step.velocity * current_step.momentum_direct_y * (t - current_step.global_time)
                z = current_step.position_z + current_step.velocity * current_step.momentum_direct_z * (t - current_step.global_time)
            else:
                x = self.track[-1].x + current_step.velocity * current_step.momentum_direct_x * TIME_STEP
                y = self.track[-1].y + current_step.velocity * current_step.momentum_direct_y * TIME_STEP
                z = self.track[-1].z + current_step.velocity * current_step.momentum_direct_z * TIME_STEP
            self.track.append(TrackPoint(t, x, y, z))
        

# ----------------------- read geant file and save geant tracks --------------------

# reads file and writes particles and their tracks to particle_list
particle_list = []    
with open(sys.argv[1], "r") as file:

    lines = text.split('\n')
    
    for line_index, line in enumerate(lines):
        
        try:
            
            if (len(line) <= 0 or line[0]=="#"):
                continue
            
            elems = line.split()
            particle = ParticleInfo()
            
            particle.picID = int(elems[0])
            particle.parentID = int(elems[1])
            particle.pic_type = elems[2]
            particle.position_x = float(elems[3])  
            particle.position_y = float(elems[4])       
            particle.position_z = float(elems[5])       
            particle.global_time = float(elems[6])      
            particle.momentum_direct_x = float(elems[7])
            particle.momentum_direct_y = float(elems[8])
            particle.momentum_direct_z = float(elems[9])
            particle.momentum_x = float(elems[10])      
            particle.momentum_y = float(elems[11])       
            particle.momentum_z = float(elems[12])       
            particle.velocity = float(elems[13])         
            particle.weight = float(elems[14])
            
            particle_list.append(particle)
            
        except:
            print("ERROR: line=\"%s\"" % (line))
            exit()


# ----------------------- compute tracks to draw using geant tracks --------------------

# shifts time of each particle (so first particle starts at time=0)
def compute_shift_time(particle_list):
    min_time = min((p.geant_steps[0].start_pi.global_time for p in particle_list))
    for p in particle_list:
        p.shift_time(min_time)
    max_time = max((p.geant_steps[-1].final_pi.global_time for p in particle_list))
    return (0., max_time)

# computes tracks to draw for each particle
def compute_tracks(particle_list, TIME_STEP, t_start, t_final):
    for p in particle_list:
        p.compute_track(TIME_STEP, t_start, t_final)


t_start, t_final = compute_shift_time(particle_list)
compute_tracks(particle_list, TIME_STEP, t_start, t_final)


# -------------- plot and animate -----------------

# borders to compute density and plot it
x_min = 0.
x_max = 25.
y_min = -2.
y_max = 2.
z_min = -10
z_max = 10.

# to plot target box
x_target_min = 10.
x_target_max = 20.
y_target_min = -5.
y_target_max = 5.
z_target_min = -5.
z_target_max = 5.

def get_density(track_index):
    N = 1000
    x_step = (x_max - x_min)/N
    z_step = (z_max - z_min)/N
    res = np.zeros(shape=(N, N))
    for p in particle_list:
        x = p.track[track_index].x
        y = p.track[track_index].y
        z = p.track[track_index].z
        if x >= x_min and x < x_max and y >= y_min and y < y_max and z >= z_min and z < z_max:
            x_index = int((x - x_min)/x_step)
            z_index = int((z - z_min)/z_step)
            res[N - z_index - 1, x_index] += 1
    return res
    
fig = plt.figure()
ax = fig.add_subplot(1,1,1)

start_im = get_density(0)
im = ax.imshow(start_im, cmap='RdBu', interpolation='none',
    extent=(x_min, x_max, z_min, z_max), animated = True)
fig.colorbar(im, ax=ax)
ax.set_xlabel("x")
ax.set_ylabel("z")

# TARGET
target_line_1, = ax.plot([x_target_min, x_target_max], [z_target_min, z_target_min], color='k', linestyle='-', linewidth=2)
target_line_2, = ax.plot([x_target_min, x_target_max], [z_target_max, z_target_max], color='k', linestyle='-', linewidth=2)
target_line_3, = ax.plot([x_target_min, x_target_min], [z_target_min, z_target_max], color='k', linestyle='-', linewidth=2)
target_line_4, = ax.plot([x_target_max, x_target_max], [z_target_min, z_target_max], color='k', linestyle='-', linewidth=2)

# animation func
def update_fig(i):
    #if i % 10 == 0:
    #    plt.savefig("pics/%d.png" % (i))
    im.set_array(get_density(i))  # get density and draw
    return im, target_line_1, target_line_2, target_line_3, target_line_4
    
ani = animation.FuncAnimation(fig, update_fig, interval=10, blit=True)

fig.tight_layout()

plt.show()


