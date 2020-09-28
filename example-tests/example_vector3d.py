import sys
sys.path.append("../bin/")
import pyHiChi as hichi

#vector3d
v = hichi.Vector3d(1.2, 2.2, 3.4)
print("volume: ", v.volume())
print("norm: ", v.norm())
print("norm2: ", v.norm2())
print("vector: ", str(v))
print(v)
print(v.x)
print(v.y)
print(v.z)
v2 = hichi.Vector3d()
print("v2: ", str(v2))
v2.x = 2.4
v2.y = 4.5
v2.z = 312
print("v2: ", str(v2))