from cmath import exp
from functools import reduce

def highest_power_dividing(l, N):
    result = 0
    while N % l == 0:
        result += 1
        N = N / l
    return result

def prime_modular_poly_mod_p(l, p):
    """Computes Phi_l in Fp, where l and p are different primes;
    
    If we assume that the order of a random elliptic curve defined over Fp is
    uniformly distributed in the Hasse interval and that the created samples are
    mostly linearly independent, we get that the complexity of the used algorithm
    is O(l^3 T) where T is the complexity of `EllipticCurve.cardinality()`."""
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
            N = E.cardinality()
            if N % l == 0:
                O = E(0, 1, 0)
                P = E.random_point()
                n = ZZ(N / l**highest_power_dividing(l, N))
                while P * n == O:
                    P = E.random_point()
                P = P * n
                while P * l != O:
                    P = P * l
                f = E.isogeny(P)
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
                print("found samples: " + str(len(samples)))
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
                                result += x**e * y**f * P(ker_gen[i])
                            else:
                                result += x**e * y**f * P(ker_gen[i]) + y**e * x**f * P(ker_gen[i])
                            i += 1
                    return result / result.coefficient({ x: degree, y: 0 })
                else:
                    target_sample_count = Integer((target_sample_count * 3/2).floor())