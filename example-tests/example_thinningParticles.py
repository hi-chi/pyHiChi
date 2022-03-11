import sys
sys.path.append("../bin/")
import pyHiChi as hichi
from random import random

def new_random_array(size) :
    p_array = hichi.ParticleArray()  #type = ELECTRON
    for i in range(size) :
        pos = hichi.Vector3d(random(),random(), random())
        mo = hichi.Vector3d(random(),random(), random())
        new_p = hichi.Particle(pos, mo, random(), hichi.ELECTRON)
        p_array.add(new_p)
    return p_array

thinning = hichi.Thinout()

#simple thinning
p_array = new_random_array(1000)
print(p_array.size())
thinning.simple_thinning(p_array, 300)
print(p_array.size())

#leveling thinning
p_array = new_random_array(1000)
print(p_array.size())
thinning.leveling_thinning(p_array)
print(p_array.size())

#number conservative thinning
p_array = new_random_array(1000)
print(p_array.size())
thinning.number_conservative_thinning(p_array, 300)
print(p_array.size())

#energy conservative thinning
p_array = new_random_array(1000)
print(p_array.size())
thinning.energy_conservative_thinning(p_array, 300)
print(p_array.size())

#energy conservative thinning
p_array = new_random_array(1000)
print(p_array.size())
features = set([hichi.Conserve.momentum,  hichi.Conserve.position,
    hichi.Conserve.energy, hichi.Conserve.dis_momentum, 
    hichi.Conserve.dis_position, hichi.Conserve.dis_energy])
thinning.conservative_thinning(p_array, 300, features)
print(p_array.size())

