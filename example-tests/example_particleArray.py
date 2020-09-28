import sys
sys.path.append("../bin/")
import pyHiChi as hichi


#ParticleArray - legacy, use Ensemble
p_array = hichi.ParticleArray()  #type = Electron
for i in range(11) :
    pos = hichi.Vector3d(1.2*i, 3.4*i, 5.6*i)
    mo = hichi.Vector3d(9.8*i, 7.6*i, 54.3*i)
    new_p = hichi.Particle(pos, mo, 0.5, hichi.ELECTRON)
    p_array.add(new_p)
print(p_array.size())
print(p_array.get_type())
print(p_array[5].get_position())
p_array.delete(5)
print(p_array[5].get_position())

for el in p_array :
    print(el.get_position())
print(p_array.size())
