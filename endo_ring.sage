from math import floor, ceil
from queue import PriorityQueue

from numpy import mintypecode

# there is a theorem that says that for a supersingular Elliptic Curve E
# defined over Fq, q = p^2 such that the (p^r)-th power Frobenius endomorphism has
# trace +/- 2p^r, we have E(F(q^r)) ~ Z/(p^r -/+ 1)Z x Z/(p^r -/+ 1)Z.
# Thus E(F(q^r)) = E[p -/+ 1]
# E.g. DeFeo, Thm 54

# it follows that if the q-th power Frobenius endomorphism has trace +/- 2p,
# then E(F(q)) ~ Z/(p^r -/+ 1)Z x Z/(p^r -/+ 1)Z and for r > 1 have 
# E(F(q^r)) ~ Z/(p^r - 1)Z x Z/(p^r - 1)Z

# We call these curves nice

def is_nice(E):
    """Returns true if E is a supersingular Elliptic Curve defined over Fq, q = p^2
    such that the q-th power Frobenius endomorphism has trace +/- 2p"""
    if not E.is_supersingular():
        return False
    p = E.base_field().characteristic()
    if not E.base_field() == GF(p**2):
        return False
    n = E.cardinality()
    return n == (p + 1)**2 or n == (p - 1)**2

def frob_trace(E):
    """Computes the trace of the q-th power frobenius endomorphism for a nice Elliptic Curve E"""
    p = E.base_field().characteristic()
    return p**2 + 1 - E.cardinality()

def rational_point_count(E, r):
    """Computes the number #E(F(q^r)), assuming that E is a nice Elliptic Curve"""
    if r == 0:
        return 1
    tr = frob_trace(E)
    q = E.base_field().characteristic()**2
    alpha = tr/2 + sqrt(CC(tr**2/4 - q))
    beta = tr/2 - sqrt(CC(tr**2/4 - q))
    return Integer((q**r + 1 - alpha**r - beta**r).real())

def rational_point_gens(E, r):
    """Computes generators of the group E(F(q^r)) for a nice Elliptic Curve E"""
    q = E.base_field().characteristic()**2
    return E.base_extend(GF(q**r)).gens()

def torsion_group_gens(E, l, rational_gens = [0, None, None]):
    """Computes generators of the group E[l] for a nice Elliptic Curve E. This uses
    the generators of the group E(F(q^r)) for different r. In case some of this data
    is already available, you can pass `rational_gens = [r, P, Q]` where P and Q are
    F(q^r)-rational points that generate E(F(q^r)). If new [r, P, Q] are computed
    during the function call, these are again stored in `rational_gens`"""
    r, P, Q = rational_gens
    if rational_point_count(E, r) % l**2 != 0:
        r += 1
        while rational_point_count(E, r) % l**2 != 0:
            r += 1
        P, Q = rational_point_gens(E, r)
    rational_gens[0] = r
    rational_gens[1] = P
    rational_gens[2] = Q

    n = sqrt(rational_point_count(E, r))

    return (P * Integer(n / l), Q * Integer(n / l))

def two_dim_discrete_log_baby_giant_step(P, Q1, Q2, n):
    """Assumes that Q1, Q2 are a basis of a free Z/nZ-module containing P.
    Computes a, b such that P = a Q1 + b Q2. This uses the baby-step-giant-step
    approach, so if n is not a prime, it will be faster to do it on each prime
    factor instead."""
    E = P.curve()
    assert P * n == E(0, 1, 0)
    assert Q1 * n == E(0, 1, 0)
    assert Q2 * n == E(0, 1, 0)
    steps = int(sqrt(n)) + 1
    R1 = Q1 * steps
    R2 = Q2 * steps
    giant_steps = { R1 * i + R2 * j: (i, j) for j in range(steps) for i in range(steps) }
    for k in range(steps):
        for l in range(steps):
            R = P + Q1 * k + Q2 * l
            if R in giant_steps:
                i, j = giant_steps[R]
                return (i * steps - k, j * steps - l)

def two_dim_discrete_log_prime_power(P, Q1, Q2, p, power):
    """Assumes that Q1, Q2 are a basis of a free Z/nZ-module containing P.
    Computes a, b such that P = a Q1 + b Q2. This assumes that n = p^power
    is a prime power, and will use power applications of the baby-step-giant-step
    method"""
    assert power >= 1
    if power == 1:
        return two_dim_discrete_log_baby_giant_step(P, Q1, Q2, p)

    # the basic idea is the exact sequence
    #   0 -> Z/pZ -> Z/p^(r+1)Z -> Z/p^rZ -> 0
    # where the two inner maps are multiplication by p^r and
    # mod p^r, respectively
    right_P, right_Q1, right_Q2 = P * p, Q1 * p, Q2 * p
    right_a, right_b = two_dim_discrete_log_prime_power(right_P, right_Q1, right_Q2, p, power - 1)
    left_P = P - right_a * Q1 - right_b * Q2
    # left_P is in the kernel of Z/p^(r+1)Z -> Z/p^rZ, thus in the image of Z/pZ -> Z/p^(r+1)Z
    left_Q1 = Q1 * p**(power - 1)
    left_Q2 = Q2 * p**(power - 1)
    left_a, left_b = two_dim_discrete_log_baby_giant_step(left_P, left_Q1, left_Q2, p)
    return (left_a * p**(power - 1) + right_a, left_b * p**(power - 1) + right_b)

def two_dim_discrete_log(P, Q1, Q2, n):
    """Assumes that Q1, Q2 are a basis of a free Z/nZ-module containing P.
    Computes a, b such that P = a Q1 + b Q2"""
    a, b = (0, 0)
    m = 1
    for (p, power) in factor(n):
        k = p**power
        a1, b1 = two_dim_discrete_log_prime_power(P * Integer(n/k), Q1 * Integer(n/k), Q2 * Integer(n/k), p, power)
        a = crt(a, a1, m, k)
        b = crt(b, b1, m, k)
        m *= k
    assert m == n
    return (a, b)

def find_nontrivial_endo_cycle(E, l):
    """Finds a nontrivial cycle in the l-isogeny graph of E and returns a list
    of l-isogenies and isomorphisms of this circle"""
    F = E.base_field()
    assert l != F.characteristic()
    open = PriorityQueue()
    found = {}
    open.put((0, E))
    found[E.j_invariant()] = None
    while True:
        prio, curve = open.get()
        for f in curve.isogenies_prime_degree(l):
            target = f.codomain().j_invariant()
            if target in found:
                if f.dual() == found[curve.j_invariant()]:
                    continue

                f = f.domain().isogeny(f.kernel_polynomial(), degree = f.degree(), codomain = found[target].codomain())
                circle_len = 1
                result = [f.dual()]
                current = curve.j_invariant()
                while found[current] != None:
                    result.append(found[current].dual())
                    circle_len += 1
                    current = found[current].domain().j_invariant()
                result.reverse()
                current = target
                while found[current] != None:
                    result.append(found[current])
                    circle_len += 1
                    current = found[current].domain().j_invariant()
                return (result, l**circle_len)
            else:
                found[target] = f
                open.put((prio + 1, target))

def eval_endo_cycle(endo_cycle, P):
    F = P.curve().base_field()
    for f in reversed(endo_cycle):
        P = f.domain().base_extend(F).isogeny(f.kernel_polynomial(), degree = f.degree(), codomain = f.codomain().base_extend(F))(P)
    return P

def min_repr(a):
    n = a.parent().order()
    a = ZZ(a)
    if a > n/2:
        a -= n
    return a

def find_mod_n_transform(E, endo, r):
    """We know that E(F(q^r)) ~ Z/nZ x Z/nZ, so an endomorphism
    of E gives a 2x2 matrix over Z/nZ that describes its action
    on n-torsion points. This function finds that matrix, by only
    evaluating the endomorphism on some points."""
    n = Integer(sqrt(rational_point_count(E, r)))
    Q1, Q2 = rational_point_gens(E, r)
    P1 = endo(Q1)
    P2 = endo(Q2)
    r1, s1 = two_dim_discrete_log(P1, Q1, Q2, n)
    assert P1 == r1 * Q1 + s1 * Q2
    r2, s2 = two_dim_discrete_log(P2, Q1, Q2, n)
    assert P2 == r2 * Q1 + s2 * Q2

    R = Integers(n)
    A = Matrix(R, [[r1, r2], [s1, s2]])
    return A

F = GF(37**2)
E = EllipticCurve(F, j = 3 + sqrt(F(15)))
assert E.is_supersingular()

I = Matrix(ZZ, [[1, 0], [0, 1]])

cycle, deg = find_nontrivial_endo_cycle(E, 3)
endo = lambda P: eval_endo_cycle(cycle, P)
A = find_mod_n_transform(E, endo, 2)
print(A - I)
print((A - I).charpoly())
print()

cycle, deg = find_nontrivial_endo_cycle(E, 5)
endo = lambda P: eval_endo_cycle(cycle, P)
B = find_mod_n_transform(E, endo, 2)
print(B + I)
print((B + I).charpoly())
print()

cycle, deg = find_nontrivial_endo_cycle(E, 7)
endo = lambda P: eval_endo_cycle(cycle, P)
C = find_mod_n_transform(E, endo, 2)
print(C - 5 * I)
print((C - 5 * I).charpoly())
print()
