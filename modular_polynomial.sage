from functools import reduce

def highest_power_dividing(l, N):
    result = 0
    while N % l == 0:
        result += 1
        N = N / l
    return result

db = ClassicalModularPolynomialDatabase()

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

    phi = db[l].change_ring(F)
    j1, j2 = phi.parent().gens()
    phi_vec = vector(F, ZZ((degree + 1) * (degree + 2) / 2))
    i = 0
    for e in range(degree + 1):
        for f in range(e + 1):
            phi_vec[i] = phi.coefficient({ j1: e, j2: f })
            i += 1
            
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

def prime_power_modular_polynomial_mod_p(l, e, p):
    phi = prime_modular_poly_mod_p(l, p)
    P = phi.parent()
    P2 = PolynomialRing(P.base_ring(), ['x', 'y', 'z'])
    x, y, z = P2.gens()
    polys = [P(x - y)]
    for i in range(1, e + 1):
        f = polys[-1](x, z).resultant(phi(z, y), z)
        f = P(f)
        if i >= 2:
            while f / polys[-2] in P:
                f = P(f / polys[-2])
        polys.append(f)
    x, y = P.gens()
    return polys[-1] / polys[-1].coefficient({ x: polys[-1](x, 0).degree(), y: 0 })
