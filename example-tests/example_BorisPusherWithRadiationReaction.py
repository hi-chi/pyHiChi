import sys
sys.path.append("../bin/")
import pyHiChi as hichi
import numpy as np

def value_E_analytical(pos, t):
    E = hichi.Vector3d(10**-5, 10**-5, 10**-5) #sin(pos.x)
    return E
    
def value_B_analytical(pos, t):
    B = hichi.Vector3d(0, 0, 0)
    return B
    
t = 0
p_array = hichi.ParticleArray()
fields_array = []
for i in range(11) :
    pos = hichi.Vector3d(1.2*i, 3.4*i, 5.6*i)
    mo = hichi.Vector3d(i*10, 0, 0)
    new_p = hichi.Particle(pos, mo, 0.5, hichi.ELECTRON)
    p_array.add(new_p)
    fields_array.append(hichi.FieldValue(value_E_analytical(pos, t), value_B_analytical(pos, t)))
    
#Boris Pusher with RadiationReaction
dt = 0.001
pusher = hichi.BorisPusher()
rr = hichi.RadiationReaction()
for k in range(11) :
    print(p_array[0].get_momentum())
    pusher(p_array, fields_array, dt)
    rr(p_array, fields_array, dt)
    t = dt * k
    for j in range(11) :
        fields_array[i].set_E(value_E_analytical(pos, t))
        fields_array[i].set_B(value_B_analytical(pos, t))

