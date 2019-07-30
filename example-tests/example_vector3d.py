import pyHiChi as pfc

#vector3d
v = pfc.vector3d(1.2, 2.2, 3.4)
print("volume: ", v.volume())
print("norm: ", v.norm())
print("norm2: ", v.norm2())
print("vector: ", v.toString())
print(v)
print(v.x)
print(v.y)
print(v.z)
v2 = pfc.vector3d()
print("v2: ", v2.toString())
v2.x = 2.4
v2.y = 4.5
v2.z = 312
print("v2: ", v2.toString())