# a small number of photons move from center

import sys
sys.path.append("../../bin/")
import pyHiChi as hichi
import numpy as np
import array

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
num_particles = 20  # 100000
particle_factor = 1
density = num_particles*particle_factor/photon_volume

electron_rest_energy = hichi.ELECTRON_MASS*hichi.c*hichi.c
average_photon_energy = 20*electron_rest_energy
energy_range = 20*electron_rest_energy
average_photon_momentum = 20*hichi.ELECTRON_MASS*hichi.c

# model parameters
time_step = wavelength/hichi.c/4.0
n_iter = 15

# particle generating

def particle_density(x, y, z):
    return density * block(np.sqrt(x*x + y*y + z*z), 0.0, photon_R)

def generate_particle(pos):
    x, y, z = (pos.x, pos.y, pos.z)
    momentum = pos / np.sqrt(x*x+y*y+z*z)*(average_photon_energy + energy_range/photon_R*np.sqrt(z*z))/hichi.c
    factor = particle_factor
    return hichi.Particle(pos, momentum, factor, hichi.ParticleTypes.PHOTON)

def generate_particles(particles):
    # walk trought all grid cells
    for i in range(grid_size_int[0]):
        for j in range(grid_size_int[1]):
            for k in range(grid_size_int[2]):
                
                # compute correct particle number in the cell using particle density
                dv_start = min_coords + hichi.Vector3d(i, j, k)*grid_step
                dv_end = dv_start + grid_step
                dv = dv_end - dv_start
                center = (dv_start + dv_end) / 2.0
                expectedParticleNum = particle_density(center.x, center.y, center.z) * dv.volume() / particle_factor
                
                # round expectedParticleNum to int
                particleNum = int(expectedParticleNum)
                if np.random.rand(1) < expectedParticleNum - float(particleNum):
                    particleNum += 1
                
                # generate the computed particle number inside the cell with random positions
                rn = np.random.rand(particleNum, 3)
                positions = [hichi.Vector3d(rn[pi][0], rn[pi][1], rn[pi][2])*dv + dv_start for pi in range(particleNum)]
                for pi in range(particleNum):
                    particles.add(
                        generate_particle(positions[pi])
                    )
    
    return particles

particles = hichi.ParticleArray(hichi.ParticleTypes.PHOTON)
absorbed_particles = hichi.ParticleArray(hichi.ParticleTypes.PHOTON)
absorbed_times = []

generate_particles(particles)

# particle absorbing, when they move out of the box

def absorbe_particles(iter):
    time = iter*time_step
    
    # filter absorbed particles
    cond = lambda ip: not min_coords <= ip[1].get_position() < max_coords
    curr_absorbed_indices = [index for index, p in filter(cond, enumerate(particles))]
    
    # save absorbed particles to the 'absorbed_particle' array
    for i in curr_absorbed_indices:
        absorbed_particles.add(particles[i])
    absorbed_times.extend([time for i in curr_absorbed_indices])
    
    # delete absorbed particles from the 'particles' array
    curr_absorbed_indices.sort()   
    for index in curr_absorbed_indices[::-1]:
        particles.delete(index)  # delete = swap with the last and size-=1

# simulation start

field = hichi.PSATDField(grid_size, min_coords, grid_step, time_step)  # field is always zero
pusher = hichi.BorisPusher()

def update(iter):
    # absorbing
    absorbe_particles(iter)
    
    # pusher
    fields_array = [hichi.FieldValue(field.get_E(p.get_position()), field.get_B(p.get_position())) for p in particles]
    pusher(particles, fields_array, time_step)
    
    # solver
    field.update_fields()
    
    # TODO: current deposition
    
for iter in range(n_iter):
    update(iter)

# print absorbed particles to file in a correct format
    
def write_geant_output(file_name):
    with open(file_name, 'wb') as file:
        file.write(bytes(particle_types[hichi.ParticleTypes.PHOTON] + "\n", encoding="utf-8"))
        file.write(bytes(str(absorbed_particles.size()) + "\n", encoding="utf-8"))
        for p, t in zip(absorbed_particles, absorbed_times):
            list_values = [
                p.get_mass(), p.get_charge(), p.get_weight(),
                p.get_position().x, p.get_position().y, p.get_position().z,
                p.get_velocity().x, p.get_velocity().y, p.get_velocity().z,
                p.get_momentum().x, p.get_momentum().y, p.get_momentum().z,
                t
            ]
            array_values = array.array("d", list_values)
            array_values.tofile(file)
            
write_geant_output(FILE_NAME)
