from itertools import count


p = 37
F = GF(p**2)
z = F.gen()
E = EllipticCurve(F, j = 1 + 2*z)
E2 = EllipticCurve(F, j = (1 + 2*z)**p)