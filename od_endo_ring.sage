from re import X
from sage.groups.generic import discrete_log
from sage.arith.misc import CRT
from discrete_log_sage import two_dim_discrete_log

def division_polynomials(E, n):
    P.<x> = E.base_field()[]
    if n == 1:
        return (x, f(x), 1)

    assert E.a1() == 0
    assert E.a3() == 0
    assert E.a2() == 0
    A = E.a4()
    B = E.a6()
    x = x
    f = x**3 + A * x + B

    x1 = 18*x**4 - 16*x*f(x) + 12*x**2*A + 2*A**2
    y1_y = -27*x**6 - 27*A*x**4 + 28*x**3*f(x) - 9*A**2*x**2 + 4*A*x*f(x) - A**3 - 8*B*f(x)
    z1 = 8*f(x)
    assert f(x1/z1) - y1_y**2/f(x)/z1**2 == 0

    for _ in range(n - 2):
        y1_squared = P(f(x1/z1)*z1**3)
        x3 = (y1_squared - 2 * y1_y * z1**2 + z1**3 * f(x)) - (x1 - z1 * x)**2 * (z1 * x + x1)
        z3 = (x1 - z1 * x)**2 * z1
        (x1, y1_y, z1) = (
            x3 * (x1 - z1 * x),
            ((z1 * f(x) - y1_y) * x3) + (z3 * (y1_y * x - f(x) * x1)),
            z3 * (x1 - z1 * x)
        )
        assert x1 in P and y1_y in P and z1 in P
        assert f(x1/z1) - y1_y**2/f(x)/z1**2 == 0
    d = gcd(x1, gcd(y1_y, z1))
    return (P(x1/d), P(y1_y/d), P(z1/d))

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
    O = E(0, 1, 0)
    gens = E.gens()
    if len(gens) == 2:
        P, Q = gens
    else:
        P = gens[0]
        Q = O
    N = E.cardinality()
    n = max(i for i in range(highest_power_dividing(l, N) + 1) if Integer(N/l**i) * P == O)
    m = max(i for i in range(highest_power_dividing(l, N) + 1) if Integer(N/l**i) * Q == O)
    assert n + m == highest_power_dividing(l, N)
    existing_l_e_torsion_points = l**(min(n, e) + min(m, e))
    # I strongly suspect that this degree is enough, even though I have not proven it yet
    degree = l**(2 * e) - existing_l_e_torsion_points

    Eext = E.base_extend(GF(q**degree))
    O = Eext(0, 1, 0)
    N = E.cardinality(extension_degree = degree)

    def random_order_le_point():
        n = coprime_part(l, N)
        P = Eext.random_point() * n
        while P * l**(e - 1) == O:
            P = Eext.random_point() * n
        while P * l**e != O:
            P = P * l
        return P

    while True:
        P = random_order_le_point()
        Q = random_order_le_point()
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
    # the discriminant of the frobenius order
    D = trace**2 - 4 * q
    # the discriminant of the maximal order
    d = D.squarefree_part()
    if (d - 1) % 4 != 0:
        d *= 4
    frobenius_order_conductor = sqrt(D/d)

    conductor = 1
    # it is easier to work with (2*frob - trace) than with (frob - trace/2),
    # so have an additional factor of 2
    for (p, power) in factor(frobenius_order_conductor * 2):
        P, Q = torsion_group_gens(E, p, power)
        for i in range(1, power + 1):
            if 2 * frob(P) == trace * P and 2 * frob(Q) == trace * Q:
                break
            conductor *= p
            P = P * p
            Q = Q * p
    K = E.frobenius_order().number_field()
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
print(division_polynomials(E, 3))
print(endo_ring(E).gens())