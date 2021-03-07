import sys
sys.path.append("../bin/")
import pyHiChi as hichi

E = hichi.Vector3d(1.2, 2.2, 3.4)
B = hichi.Vector3d(34, 5.6, 7.8)

#FieldValue
f = hichi.FieldValue(E, B)
print(f.get_E())
print(f.get_B())
f.set_E(B)
f.set_B(E)
print(f.get_E())
print(f.get_B())

f2 = hichi.FieldValue(E.z, E.y, E.x, B.x, B.x, B.x) #Ex, Ey, Ez, Bx, By, Bz
print(f2.get_E())
print(f2.get_B())