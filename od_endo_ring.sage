from re import X
from sage.groups.generic import discrete_log
from sage.arith.misc import CRT
from discrete_log_sage import two_dim_discrete_log
from division_polynomial_sage import division_polynomials

def highest_power_dividing(l, N):
    result = 0
    while N % l == 0:
        result += 1
        N = N / l
    return result

def coprime_part(l, N):
    return Integer(N / l**highest_power_dividing(l, N))

def torsion_group_gens(E, l, e):
    """Computes generators of E[l^e], where l is a prime"""
    if "torsion_group_gens_cache" not in E:
        E.torsion_group_gens_cache = {}
    elif l**e in E.torsion_group_gens_cache:
        return E.torsion_group_gens_cache[l**e]

    assert l.is_prime()
    q = E.base_field().order()
    
    f, _, h = division_polynomials(E, l**e)
    P = PolynomialRing(E.base_field(), ['x', 'z'])
    x, z = P.gens()
    h = P(h / gcd(h, f))
    h = P(P(h)(x = x/z) * z**h.degree()) * z
    degree = 1
    for (c, _) in factor(h):
        degree = lcm(degree, c.degree())
    # we need to take the root of x^3 + Ax + B
    degree = degree * 2

    F = GF(q**degree)
    Eext = E.base_extend(F)
    O = Eext(0, 1, 0)

    def points():
        for (u, _) in factor(h):
            if u == z:
                yield O
                continue
            for (v, _) in factor(u.change_ring(F)):
                t, _ = v.coefficients()
                t = -1/t
                r = sqrt(t**3 + E.a4()*t + E.a6())
                yield Eext(t, r)
    
    point_iter = points()
    P = next(point_iter)
    if P == O:
        P = next(point_iter)
    assert P != O and P * l**e == O
    while True:
        Q = next(point_iter)
        try:
            discrete_log(P * l**(e - 1), Q, ord = l**e, operation = '+')
        except ValueError:
            E.torsion_group_gens_cache[l**e] = (P, Q)
            return (P, Q)

def endo_ring(E):
    assert E.is_ordinary()
    # first, we find an isomorphic elliptic curve over as small a base field as possible
    j = E.j_invariant()
    p = E.base_field().characteristic()
    degree = j.minpoly().degree()
    q = p**degree
    F = GF(q)
    E = EllipticCurve(F, j = j)
    frob = lambda P: P.curve()(P[0]**q, P[1]**q, P[2]**q)

    trace = E.trace_of_frobenius()
    K = E.frobenius_order().number_field()
    # it is easier to work with phi := (2*frob - trace) than with frob or frob - trace/2
    phi_order_conductor = K.order(E.frobenius() * 2 - trace).index_in(K.maximal_order())

    conductor = 1
    for (p, power) in factor(phi_order_conductor):
        P, Q = torsion_group_gens(E, p, power)
        for i in range(1, power + 1):
            if 2 * frob(P) == trace * P and 2 * frob(Q) == trace * Q:
                break
            conductor *= p
            P = P * p
            Q = Q * p
    alpha, = K.maximal_order().ring_generators()
    return K.order(alpha * conductor)

def kernel(A):
    """Computes the kernel of a matrix A defined over Z/mZ"""
    R = A.base_ring()
    m = R.characteristic()
    A = Matrix(ZZ, [[ZZ(a) for a in col] for col in A])
    B, U, V = A.smith_form()
    n = len(B.rows())
    result = Matrix.identity(R, n)
    for i in range(n):
        d = B[i, i]
        result[i] = V * result[i] * Integer(lcm(m, d)/d)
    result = result.transpose()
    return result

def prime_powers_leq(n):
    for l in Primes():
        if l <= n:
            yield (l, Integer(log(n, l)))
        else:
            return

class Endo:

    def __init__(self, j, a, b):
        """The endomorphism a + bf where f is the q-th power frobenius endomorphism"""
        self.p = j.parent().characteristic()
        self.q = p**j.minpoly().degree()
        self.F = GF(self.q)
        self.E = EllipticCurve(self.F, j = j)
        self.a = a
        self.b = b

    def deg(self):
        # this follows from the fact that #E(Fq) = deg(1 - pi)
        return self.a**2 + self.b**2 * self.q + 2 * self.a * self.b * (self.q + 1 - self.E.cardinality())

    def get_prime_power_transform(self, l, e):
        if l**e in self.transforms:
            return self.transforms[l**e]
        P, Q = torsion_group_gens(self.E, l, e)
        a, b = two_dim_discrete_log(self.eval(P), P, Q, l**e)
        c, d = two_dim_discrete_log(self.eval(Q), P, Q, l**e)
        R = Integers(l**e)
        self.transforms[l**e] = Matrix(R, [[a, c], [b, d]])
        return self.transforms[l**e]

    def kernel(self):
        P0, Q0 = self.E(0, 1, 0), self.E(0, 1, 0)
        n0, m0 = 1, 1
        N = 1
        for (l, e) in factor(self.deg()):
            P, Q = torsion_group_gens(self.E, l, e)
            transform = self.get_prime_power_transform(l, e) + rhs.get_prime_power_transform(l, e)
            A = kernel(transform)
            P1 = ZZ(A[0, 0]) * P + ZZ(A[0, 1]) * Q
            Q1 = ZZ(A[1, 0]) * P + ZZ(A[1, 1]) * Q
            n1 = gcd(ZZ(A[0, 0]), ZZ(A[0, 1]))
            m1 = gcd(ZZ(A[1, 0]), ZZ(A[1, 1]))
            n0, s, t = xgcd(n0 * l**e, n1 * N)
            P0 = s * P0 + t * P1
            m0, s, t = xgcd(m0 * l**e, m1 * N)
            Q0 = s * Q0 + t * Q1
            N = N * l**e
        return (P0, Q0)

        
F = GF(37**2)
E = EllipticCurve(F, j = 3 * F.gen())
print(endo_ring(E).gens())