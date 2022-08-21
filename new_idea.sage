from modular_polynomial_sage import prime_modular_poly_mod_p

p = next_prime(12)
polys = []
poly_ring = PolynomialRing(GF(p), 'x')
for l in [3, 5]:
    phi = prime_modular_poly_mod_p(l, p)
    phi2 = prime_power_modular_poly_mod_p(l, 2, p)
    P = PolynomialRing(GF(p), ['x', 'y1', 'y2', 'y3', *['b' + str(i) for i in range(l + 1)]])
    x = P.gen(0)
    y1, y2, y3 = P.gen(1), P.gen(2), P.gen(3)
    bs = P.gens()[(l + 2):]
    I = P.ideal([
        phi(x, y1),
        phi(x, y2),
        phi(x, y3),
        phi2(y1, y2),
        phi2(y1, y3),
        phi2(y2, y3),
        sum(y1**i * bs[i] for i in range(l + 1)) - 1,
        sum(y2**i * bs[i] for i in range(l + 1)) - 1,
        sum(y3**i * bs[i] for i in range(l + 1)) - 1
    ])
    J = I.elim([y1, y2, y3])

    a = phi.polynomial(phi.parent().gen(1))
    A = Matrix(poly_ring, l + 1, l + 1)
    for i in range(l + 1):
        A[i, 0] = a[i]
    for i in range(l):
        A[i + 1, i] = 1
    b = A.pow(p^2 - 1 - l)[:, 0]

    for f in J.gens():
        polys.append(f({ x: poly_ring.gen(), **{ bs[i]: b[i] for i in range(l + 1) } }))

for p in polys:
    print(p)