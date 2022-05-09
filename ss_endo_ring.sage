from isogeny_graph_sage import isogeny_graph_bfs
from discrete_log_sage import two_dim_discrete_log

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
    alpha = tr/2 + sqrt(tr**2/4 - q)
    beta = tr/2 - sqrt(tr**2/4 - q)
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


def find_nontrivial_endo_cycle(E, l):
    """Finds a nontrivial cycle in the l-isogeny graph of E and returns a list
    of l-isogenies and isomorphisms of this circle"""
    cycles = isogeny_graph_bfs(E, isogeny_degree = l, find_cycle = True)
    return next(cycles), next(cycles), next(cycles)

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

    A = Matrix(QQ, [[r1/n, r2/n], [s1/n, s2/n]])
    return A, n

F = GF(37**2)
E = EllipticCurve(F, j = 3 + sqrt(F(15)))
assert E.is_supersingular()

I = Matrix(ZZ, [[1, 0], [0, 1]])

cycle, cycle2, cycle3 = find_nontrivial_endo_cycle(E, 3)
endo = lambda P: eval_endo_cycle(cycle, P)
A, n = find_mod_n_transform(E, endo, 2)
print(A)
A, n = find_mod_n_transform(E, endo, 6)
print(A)
print()