import pythonModule as pfc
import random

def newRandomArray(size) :
	pArray = pfc.particleArray()  #type = Electron
	for i in range(size) :
		pos = pfc.vector3d(1.2*i, 3.4*i, 5.6*i)
		mo = pfc.vector3d(9.8*i, 7.6*i, 54.3*i)
		newP = pfc.particle(pos, mo, 0.5, pfc.Electron)
		pArray.add(newP) #or pArray.pushBack
	return pArray
	
#simple thinning
pArray =  newRandomArray(1000)
print(pArray.size())
pfc.simpleThinning(pArray, 300)
print(pArray.size())

#leveling thinning
pArray =  newRandomArray(1000)
print(pArray.size())
pfc.levelingThinning(pArray)
print(pArray.size())

#number conservative thinning
pArray =  newRandomArray(1000)
print(pArray.size())
pfc.numberConservativeThinning(pArray, 300)
print(pArray.size())

#energy conservative thinning
pArray =  newRandomArray(1000)
print(pArray.size())
pfc.energyConservativeThinning(pArray, 300)
print(pArray.size())

