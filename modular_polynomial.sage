from cmath import exp
from functools import reduce

def prime_modular_poly_mod_p(l, p):
    """Computes Phi_l in Fp, where l and p are different primes"""
    assert l != p
    assert l.is_prime()
    assert p.is_prime()
    degree = l + 1
    F = GF(p)
    freedom = (degree + 1) * (degree + 2) / 2
    # we include symmetric samples
    target_sample_count = (degree + 1)**2
    samples = []
    for ext_deg in range(1, 6):
        F = GF(p**ext_deg)
        for j1 in F:
            if j1.minpoly().degree() != ext_deg:
                continue
            E = EllipticCurve(F, j = j1)
            n = E.cardinality()
            if n % l == 0:
                O = E(0, 1, 0)
                P = E.random_point()
                while P * ZZ(n/l) == O:
                    P = E.random_point()
                f = E.isogeny(P * ZZ(P.order()/l))
                j2 = f.codomain().j_invariant()

                sample = vector(F, ZZ(freedom))
                i = 0
                for e in range(degree + 1):
                    for f in range(e + 1):
                        if e == f:
                            sample[i] = j1**e * j2**f
                        else:
                            sample[i] = j1**e * j2**f + j1**f * j2**e
                        i += 1

                samples.append(sample)

            if len(samples) >= target_sample_count:
                K = GF(p**reduce(lcm, range(1, ext_deg + 1)))
                A = Matrix(K, samples)
                ker = kernel(A.transpose())
                if ker.dimension() == 1:
                    ker_gen, = ker.basis()
                    result = 0
                    P = PolynomialRing(GF(p), ['x', 'y'])
                    x, y = P.gens()
                    i = 0
                    for e in range(degree + 1):
                        for f in range(e + 1):
                            if e == f:
                                result += x**e * y**f * ker_gen[i]
                            else:
                                result += x**e * y**f * ker_gen[i] + y**e * x**f * ker_gen[i]
                            i += 1
                    return result / result.coefficient({ x: degree, y: 0 })
                else:
                    target_sample_count = Integer((target_sample_count * 3/2).floor())

p = 37
l = 7
x, y = PolynomialRing(GF(p), ['x', 'y']).gens()
actual = prime_modular_poly_mod_p(l, p)
print(actual == ClassicalModularPolynomialDatabase()[l](x, y))