import sys
sys.path.append("../bin/")
import pyHiChi as hichi


# ensemble
ensemble = hichi.Ensemble()
for i in range(11) :
    pos = hichi.Vector3d(1.2*i, 1.3*i, 1.6*i)
    mo = hichi.Vector3d(1.1*i, 1.4*i, 1.5*i)
    new_p = hichi.Particle(pos, mo, 0.5, hichi.ELECTRON)
    ensemble.add(new_p)
for i in range(10) :
    pos = hichi.Vector3d(13*i, 14*i, 17*i)
    mo = hichi.Vector3d(12*i, 15*i, 16*i)
    new_p = hichi.Particle(pos, mo, 0.5, hichi.POSITRON)
    ensemble.add(new_p)
for i in range(13) :
    pos = hichi.Vector3d(140*i, 150*i, 180*i)
    mo = hichi.Vector3d(130*i, 160*i, 170*i)
    new_p = hichi.Particle(pos, mo, 0.5, hichi.PROTON)
    ensemble.add(new_p)

print('Count Particles: ', ensemble.size())
print('Count Electron: ', ensemble['Electron'].size()) #use index 'Electron' or hichi.ELECTRON
print('Count Positron: ', ensemble['Positron'].size()) #same Electron
print('Count Proton: ', ensemble['Proton'].size())     #same Electron
    

print('Positions Electron: ')
for elem in ensemble[hichi.ELECTRON] :
    print(elem.get_position())
positron_array = ensemble[hichi.POSITRON]
print('Position second Positron')
print(positron_array[1].get_position())
