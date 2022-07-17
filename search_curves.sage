from modular_polynomial_sage import prime_modular_poly_mod_p, prime_power_modular_polynomial_mod_p

p = next_prime(20)
l = 2
e = int(log(p, l)) + 1
F = GF(p**2)
P = PolynomialRing(F, ['x'])
x = P.gen()
alpha = (x**p + x).roots()[1][0]

P = PolynomialRing(F, ['x',])
x, = P.gens()

phi1 = prime_power_modular_polynomial_mod_p(l, e, p)
phi2 = prime_power_modular_polynomial_mod_p(l, e + 1 , p)
j_invariants = [a for (a, _) in gcd(phi1(x, x**p), phi2(x, x**p)).roots()]
for j in j_invariants:
    print(j, EllipticCurve(F, j = j).is_supersingular())
