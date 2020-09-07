import sys
sys.path.append("../bin/")
import pyHiChi as hichi
import random

def new_random_array(size) :
    p_array = hichi.ParticleArray()  #type = ELECTRON
    for i in range(size) :
        pos = hichi.Vector3d(1.2*i, 3.4*i, 5.6*i)
        mo = hichi.Vector3d(9.8*i, 7.6*i, 54.3*i)
        new_p = hichi.Particle(pos, mo, 0.5, hichi.ELECTRON)
        p_array.add(new_p)
    return p_array
    
#simple thinning
p_array = new_random_array(1000)
print(p_array.size())
hichi.simple_thinning(p_array, 300)
print(p_array.size())

#leveling thinning
p_array = new_random_array(1000)
print(p_array.size())
hichi.leveling_thinning(p_array)
print(p_array.size())

#number conservative thinning
p_array = new_random_array(1000)
print(p_array.size())
hichi.number_conservative_thinning(p_array, 300)
print(p_array.size())

#energy conservative thinning
p_array = new_random_array(1000)
print(p_array.size())
hichi.energy_conservative_thinning(p_array, 300)
print(p_array.size())

