import pyHiChi as pfc


#ParticleArray - legacy, use Ensemble
pArray = pfc.particleArray()  #type = Electron
for i in range(11) :
	pos = pfc.vector3d(1.2*i, 3.4*i, 5.6*i)
	mo = pfc.vector3d(9.8*i, 7.6*i, 54.3*i)
	newP = pfc.particle(pos, mo, 0.5, pfc.Electron)
	pArray.add(newP)
print(pArray.size())
print(pArray.getType())
print(pArray[5].getPosition())
pArray.delete(5)
print(pArray[5].getPosition())

for el in pArray :
    print(el.getPosition())
print(pArray.size())
