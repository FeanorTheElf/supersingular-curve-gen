from re import X
from zipfile import ZIP_BZIP2
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

    print("Finding " + str(l**e) + "-torsion generators")

    assert l.is_prime()

    def extension_field(d):
        F = E.base_field()
        P = PolynomialRing(F, ['t'])
        t = P.gen(0)
        poly = t**d
        while not poly.is_prime():
            for i in range(min(d - 1, 10)):
                poly += F.random_element() * t**i
        return F.extension(poly)

    def extension_field_sqrt(f):
        P = PolynomialRing(f.parent(), ['T'])
        T = P.gen(0)
        roots = (T**2 - f).roots()
        assert len(roots) == 2 or roots[0][1] == 2
        return roots[0][0]
    
    f, _, h = division_polynomials(E, l**e)
    P = PolynomialRing(E.base_field(), ['x'])
    h = P(h / gcd(h, f))
    degree = 1
    for (c, _) in factor(h):
        degree = lcm(degree, c.degree())

    degree = degree * 2
    print("Extension degree is " + str(degree))
    F = extension_field(degree)
    Eext = E.base_extend(F)
    O = Eext(0, 1, 0)

    def points():
        for (u, _) in factor(h):
            print("Factoring " + str(u.change_ring(F)))
            for (x, _) in u.change_ring(F).roots():
                y = extension_field_sqrt(x**3 + E.a4()*x + E.a6())
                yield Eext(x, y)
    
    point_iter = points()
    P = next(point_iter)
    while P * l**(e - 1) == O:
        P = next(point_iter)
    assert P != O and P * l**e == O
    while True:
        Q = next(point_iter)
        try:
            discrete_log(Q * l**(e - 1), P, ord = l**e, operation = '+')
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
    assert phi_order_conductor % 2 == 0

    conductor = 1
    print("Order of phi: " + str(factor(phi_order_conductor)))
    for (p, power) in factor(phi_order_conductor):
        P, Q = torsion_group_gens(E, p, power)
        for i in range(1, power + 1):
            if 2 * frob(P) == trace * P and 2 * frob(Q) == trace * Q:
                break
            conductor *= p
            P = P * p
            Q = Q * p

    if highest_power_dividing(2, conductor) == highest_power_dividing(2, phi_order_conductor):
        conductor = Integer(conductor/2)
    d = K.maximal_order().discriminant() # Z[this] is the maximal order
    return K.order((sqrt(K(d)) + d)/2 * conductor)

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
        if d != 0:
            result[i] = V * result[i] * Integer(lcm(m, d)/d)
        else:
            result[i] = V * result[i]
    result = result.transpose()
    return result

def prime_powers_leq(n):
    for l in Primes():
        if l <= n:
            yield (l, Integer(log(n, l)))
        else:
            return

class EndoRing:

    def __init__(self, j):
        """The endomorphism a + bf where f is the q-th power frobenius endomorphism"""
        self.p = j.parent().characteristic()
        self.q = self.p**j.minpoly().degree()
        self.F = GF(self.q)
        self.E = EllipticCurve(self.F, j = j)
        assert self.E.is_ordinary()
        self.fundamental_generator_isogeny = None

    def _deg(self, a, b):
        # this follows from the fact that #E(Fq) = deg(1 - pi)
        return ZZ((self.E.frobenius() * b + a).norm())

    def _generator_isogeny(self):
        trace = self.E.trace_of_frobenius()
        if trace % 2 == 0:
            phi = (-ZZ(trace/2), 1)
        else:
            phi = (-trace, 2)
        P, Q, n, m = self._kernel(*phi)
        d = gcd(n, m)
        P = P * d
        Q = Q * d
        n = ZZ(n / d)
        m = ZZ(m / d)
        _, s, t = xgcd(n, m)
        P = s * Q + m * P
        # now this is the looked for endomorphism, up to automorphism
        isogeny = self.E.isogeny(P, codomain = self.E, degree = n * m)
        # for now, do not handle the nontrivial automorphism cases
        assert self.j != 0 and self.j != 1728
        P = E.random_point()
        if isogeny(P) * d == P * phi[0] + self._eval_frob(P) * phi[1]:
            return isogeny
        else:
            return -isogeny

    def _eval_frob(self, P):
        return P.curve()(P[0]**self.q, P[1]**self.q, P[2]**self.q)

    def _get_frob_prime_power_transform(self, l, e, P, Q):
        """Expects that E[l^e] has the Z/l^eZ-module basis P, Q"""
        a, b = two_dim_discrete_log(self._eval_frob(P), P, Q, l**e)
        c, d = two_dim_discrete_log(self._eval_frob(Q), P, Q, l**e)
        R = Integers(l**e)
        return Matrix(R, [[a, c], [b, d]])

    def _kernel(self, a, b):
        P0, Q0 = self.E(0, 1, 0), self.E(0, 1, 0)
        n0, m0 = 1, 1
        N = 1
        for (l, e) in factor(self._deg(a, b)):
            P, Q = torsion_group_gens(self.E, l, e)
            transform = self._get_frob_prime_power_transform(l, e, P, Q) * b + Matrix(ZZ, [[a, 0], [0, a]])
            A = kernel(transform)
            P1 = ZZ(A[0, 0]) * P + ZZ(A[0, 1]) * Q
            Q1 = ZZ(A[1, 0]) * P + ZZ(A[1, 1]) * Q
            n1 = gcd(ZZ(A[0, 0]), ZZ(A[0, 1]))
            m1 = gcd(ZZ(A[1, 0]), ZZ(A[1, 1]))
            n0, s, t = xgcd(n0 * l**e, n1 * N)
            F = GF(self.p**lcm(P0.curve().base_field().degree(), P1.curve().base_field().degree()))
            E = self.E.base_extend(F)
            P0 = s * E(P0) + t * E(P1)
            m0, s, t = xgcd(m0 * l**e, m1 * N)
            Q0 = s * E(Q0) + t * E(Q1)
            N = N * l**e
        return (P0, Q0, N/n0, N/m0)
        
F = GF(101**2)
z = F.gen()
for j in F:
    E = EllipticCurve(F, j = j)
    if E.is_supersingular():
        continue

    K = E.frobenius_order().number_field()
    phi_order_conductor = K.order(E.frobenius() * 2 - E.trace_of_frobenius()).index_in(K.maximal_order())
    skip = False
    for (l, e) in factor(phi_order_conductor):
        if l**e > 20:
            print("Skipping j = " + str(j) + " because [O_K : Z[phi]] has large prime power factor " + str(l**e))
            skip = True
            break
    if skip:
        continue

    R = endo_ring(E)
    print(j, R.class_number())