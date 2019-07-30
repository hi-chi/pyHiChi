import pyHiChi as pfc

v = pfc.vector3d(1.2, 2.2, 3.4)
v2 = pfc.vector3d(34, 5.6, 7.8)

#Particle
p = pfc.particle()
print(p.getPosition())
p.setPosition(v)
print(p.getPosition())
p.setPosition(pfc.vector3d(1.2, 2.1, 1.2))
print(p.getPosition())
p.setMomentum(v2)
print(p.getMomentum())
p.setMomentum(pfc.vector3d(3.4, 5.6, 7.8))
print(p.getMomentum())

mom = p.getMomentum()
print(mom.x)
print(p.getMomentum().y)

p.setVelocity(v)
print(p.getVelocity())
print(p.getMomentum())

print(p.getType())
print(p.getMass())
print(p.getCharge())
print(p.getWeight())

pos = pfc.vector3d(1.2, 3.4, 5.6)
mo = pfc.vector3d(9.8, 7.6, 54.3)
newP = pfc.particle(pos, mo, 0.5, pfc.Positron) #position, momentum, weight, type 
print(newP.getPosition())
print(newP.getMomentum())
print(newP.getType())
print(newP.getMass())
print(newP.getCharge())
print(newP.getWeight())

newP2 = pfc.particle(pos, mo) #position, momentum, weight = 1, type = Electron 
print(newP2.getPosition())
print(newP2.getMomentum())
print(newP2.getType())
print(newP2.getMass())
print(newP2.getCharge())
print(newP2.getWeight())