import sys
sys.path.append("../bin/")
import pyHiChi as pfc


# Ensemble
Ensemble = pfc.ensemble()
for i in range(11) :
    pos = pfc.vector3d(1.2*i, 1.3*i, 1.6*i)
    mo = pfc.vector3d(1.1*i, 1.4*i, 1.5*i)
    newP = pfc.particle(pos, mo, 0.5, pfc.Electron)
    Ensemble.add(newP)
for i in range(10) :
    pos = pfc.vector3d(13*i, 14*i, 17*i)
    mo = pfc.vector3d(12*i, 15*i, 16*i)
    newP = pfc.particle(pos, mo, 0.5, pfc.Positron)
    Ensemble.add(newP)
for i in range(13) :
    pos = pfc.vector3d(140*i, 150*i, 180*i)
    mo = pfc.vector3d(130*i, 160*i, 170*i)
    newP = pfc.particle(pos, mo, 0.5, pfc.Proton)
    Ensemble.add(newP)

print('Count Particles: ', Ensemble.size())
print('Count Electron: ', Ensemble['Electron'].size()) #use index 'Electron' or pfc.Electron
print('Count Positron: ', Ensemble['Positron'].size()) #same Electron
print('Count Proton: ', Ensemble['Proton'].size())     #same Electron
    

print('Positions Electron: ')
for elem in Ensemble[pfc.Electron] :
    print(elem.getPosition())
positronArray = Ensemble[pfc.Positron]
print('Position second Positron')
print(positronArray[1].getPosition())