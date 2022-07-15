
def defined_fp(D, p):
    return len(gp.Qfb(1, 0, D).qfbsolve(p)) > 0

def curve(D, p):
    f = hilbert_class_polynomial(D)
    for r in range(1, 10):
        roots = f.change_ring(GF(p**r)).roots()
        if len(roots) > 0:
            E = EllipticCurve(GF(p**r), j = roots[0][0])
            return E

def is_ordinary(D, p):
    return D % p == 0

D = -4
p = 11
print(curve(D, p).frobenius_order().discriminant())