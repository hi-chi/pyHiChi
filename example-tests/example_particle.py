import sys
sys.path.append("../bin/")
import pyHiChi as hichi

v = hichi.Vector3d(1.2, 2.2, 3.4)
v2 = hichi.Vector3d(34, 5.6, 7.8)

#Particle
p = hichi.Particle()
print(p.get_position())
p.set_position(v)
print(p.get_position())
p.set_position(hichi.Vector3d(1.2, 2.1, 1.2))
print(p.get_position())
p.set_momentum(v2)
print(p.get_momentum())
p.set_momentum(hichi.Vector3d(3.4, 5.6, 7.8))
print(p.get_momentum())

mom = p.get_momentum()
print(mom.x)
print(p.get_momentum().y)

p.set_velocity(v)
print(p.get_velocity())
print(p.get_momentum())

print(p.get_type())
print(p.get_mass())
print(p.get_charge())
print(p.get_weight())

pos = hichi.Vector3d(1.2, 3.4, 5.6)
mo = hichi.Vector3d(9.8, 7.6, 54.3)
new_p = hichi.Particle(pos, mo, 0.5, hichi.POSITRON) #position, momentum, weight, type 
print(new_p.get_position())
print(new_p.get_momentum())
print(new_p.get_type())
print(new_p.get_mass())
print(new_p.get_charge())
print(new_p.get_weight())

newP2 = hichi.Particle(pos, mo) #position, momentum, weight = 1, type = Electron 
print(newP2.get_position())
print(newP2.get_momentum())
print(newP2.get_type())
print(newP2.get_mass())
print(newP2.get_charge())
print(newP2.get_weight())
