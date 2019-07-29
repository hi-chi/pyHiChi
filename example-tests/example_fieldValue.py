import pythonModule as pfc

E = pfc.vector3d(1.2, 2.2, 3.4)
B = pfc.vector3d(34, 5.6, 7.8)

#FieldValue
f = pfc.field(E, B)
print(f.getE())
print(f.getB())
f.setE(B)
f.setB(E)
print(f.getE())
print(f.getB())

f2 = pfc.field(E.z, E.y, E.x, B.x, B.x, B.x) #Ex, Ey, Ez, Bx, By, Bz
print(f2.getE())
print(f2.getB())