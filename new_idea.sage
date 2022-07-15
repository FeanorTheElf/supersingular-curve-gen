from modular_polynomial_sage import prime_modular_poly_mod_p

p = next_prime(10000)
F = GF(p**2)
P = PolynomialRing(F, ['x1', 'x2', 'x3', 'x4'])
xs = P.gens()
phis = [
    prime_modular_poly_mod_p(3, p),
    prime_modular_poly_mod_p(5, p),
    prime_modular_poly_mod_p(7, p),
    prime_modular_poly_mod_p(11, p),
    prime_modular_poly_mod_p(13, p)
]
constraints = [
    phis[0](xs[0], xs[1]),
    phis[1](xs[0], xs[2]),
    phis[2](xs[0], xs[3]),
    phis[3](xs[1], xs[3]),
    phis[4](xs[2], xs[3])
]
I = P.ideal(constraints)
points = I.variety()
print(points)