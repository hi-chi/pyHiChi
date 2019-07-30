import pyHiChi as pfc
import numpy as np

def valueEAnalytical(pos, t):
	E = pfc.vector3d(1, 0, 0) #sin(pos.x)
	return E
	
def valueBAnalytical(pos, t):
	B = pfc.vector3d(0, 0, 0)
	return B
	
t = 0
pArray = pfc.particleArray()
fieldsArray = []
for i in range(11) :
	pos = pfc.vector3d(1.2*i, 3.4*i, 5.6*i)
	mo = pfc.vector3d(i*10**16, 0, 0)
	newP = pfc.particle(pos, mo, 0.5, pfc.Electron)
	pArray.add(newP)
	fieldsArray.append(pfc.field(valueEAnalytical(pos, t), valueBAnalytical(pos, t)))
	
#Boris Pusher
dt = 0.1
pusher = pfc.BorisPusher()
for k in range(11) :
	print(pArray[1].getMomentum())
	pusher(pArray, fieldsArray, dt)
	t = dt * k
	for j in range(11) :
		fieldsArray[i].setE(valueEAnalytical(pos, t))
		fieldsArray[i].setB(valueBAnalytical(pos, t))

