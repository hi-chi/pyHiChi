# a small number of photons move from center

import sys
sys.path.append("../../bin/")
import pyHiChi as hichi
import numpy as np

FILE_NAME = 'hichi_particles.txt'

particle_types = {
    hichi.ParticleTypes.ELECTRON : "Electron",
    hichi.ParticleTypes.POSITRON : "Positron",
    hichi.ParticleTypes.PROTON : "Proton",
    hichi.ParticleTypes.PHOTON : "Photon",
}

np.random.seed(3)  # fixed seed

def sign(x):
    return -1.0 if x < 0.0 else 1.0 if x > 0.0 else 0.0

def block(x, xmin, xmax):
    return 0.5 * (sign(x - xmin) + sign(xmax - x))

eV = 1.6*10e-12
MeV = 1e6*eV

# computational grid size
grid_size_int = [4, 4, 4]
grid_size = hichi.Vector3d(grid_size_int[0], grid_size_int[1], grid_size_int[2])

# box size
wavelength = 0.9e-4

L = wavelength
size_R = 2*L

min_coords = hichi.Vector3d(-size_R, -size_R, -size_R)
max_coords = hichi.Vector3d(size_R, size_R, size_R)

grid_step = (max_coords - min_coords) / grid_size

# start photon distribution
photon_R = 1.0*wavelength

photon_volume = 4.0/3.0*np.pi*(photon_R*photon_R*photon_R)
num_particles = 20
particle_factor = 1
density = num_particles*particle_factor/photon_volume

electron_rest_energy = hichi.ELECTRON_MASS*hichi.c*hichi.c
average_photon_energy = 20*electron_rest_energy
energy_range = 20*electron_rest_energy
average_photon_momentum = 20*hichi.ELECTRON_MASS*hichi.c

# model parameters
time = 0.0

time_step = wavelength/hichi.c/4.0
n_iter = 15

# photon generating

def particle_density(x, y, z):
    return density * block(np.sqrt(x*x + y*y + z*z), 0.0, photon_R)

def generate_particle(pos):
    x, y, z = (pos.x, pos.y, pos.z)
    momentum = pos / np.sqrt(x*x+y*y+z*z)*(average_photon_energy + energy_range/photon_R*np.sqrt(z*z))/hichi.c
    factor = particle_factor
    return hichi.Particle(pos, momentum, factor, hichi.ParticleTypes.PHOTON)

def generate_particles():
    particles = hichi.ParticleArray(hichi.ParticleTypes.PHOTON)

    for i in range(grid_size_int[0]):
        for j in range(grid_size_int[1]):
            for k in range(grid_size_int[2]):
                dv_start = min_coords + hichi.Vector3d(i, j, k)*grid_step
                dv_end = dv_start + grid_step
                dv = dv_end - dv_start
                center = (dv_start + dv_end) / 2.0
                expectedParticleNum = particle_density(center.x, center.y, center.z) * dv.volume() / particle_factor
                    
                particleNum = int(expectedParticleNum)
                if np.random.rand(1) < expectedParticleNum - float(particleNum):
                    particleNum += 1
                    
                #print(dv_start, dv_end, center, expectedParticleNum)
                
                rn = np.random.rand(particleNum, 3)
                positions = [hichi.Vector3d(rn[pi][0], rn[pi][1], rn[pi][2])*dv + dv_start for pi in range(particleNum)]
                for pi in range(particleNum):
                    particles.add(
                        generate_particle(positions[pi])
                    )
    
    return particles

particles = generate_particles()

# particle absorbing, when they move out of the box

absorbed_particles = []
absorbed_times = []

def absorbe_particles(time):
    cond = lambda ip: not min_coords <= ip[1].get_position() < max_coords
    curr_absorbed_indices = [index for index, p in filter(cond, enumerate(particles))]
    
    absorbed_particles.extend([particles[i] for i in curr_absorbed_indices])
    absorbed_times.extend([time for p in curr_absorbed_indices])
    
    for index in curr_absorbed_indices[::-1]:  # curr_absorbed_indices is a sorted array
        particles.delete(index)  # delete = swap with the last and size-=1   

# simulation start

field = hichi.PSATDField(grid_size, min_coords, grid_step, time_step)  # field is always zero
pusher = hichi.BorisPusher()

def update():
    absorbe_particles(time)

    fields_array = [hichi.FieldValue(field.get_E(p.get_position()), field.get_B(p.get_position())) for p in particles]
    pusher(particles, fields_array, time_step)
    field.update_fields()
    
for iter in range(n_iter):
    update()
    time += time_step

# print absorbed particles to file in a correct format
    
def write_geant_output(file_name):
    with open(file_name, 'w') as file:
        file.write("Type Mass Charge Factor Position(x) Position(y) Position(z) " + 
                   "Velocity(x) Velocity(y) Velocity(z) Momentum(x) Momentum(y) Momentum(z) Time\n")
        for p, t in zip(absorbed_particles, absorbed_times):
            res = ""
            res += str(particle_types[p.get_type()]) + " "
            res += str(p.get_mass()) + " "
            res += str(p.get_charge()) + " "
            res += str(p.get_weight()) + " "
            res += str(p.get_position().x) + " " + str(p.get_position().y) + " " + str(p.get_position().z) + " "
            res += str(p.get_velocity().x) + " " + str(p.get_velocity().y) + " " + str(p.get_velocity().z) + " "
            res += str(p.get_momentum().x) + " " + str(p.get_momentum().y) + " " + str(p.get_momentum().z) + " "
            res += str(t) + "\n"
            file.write(res)
            
write_geant_output(FILE_NAME)
