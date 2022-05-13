from itertools import count
from queue import PriorityQueue
from sage.schemes.elliptic_curves.ell_generic import is_EllipticCurve
from random import uniform, seed

seed(int(0))

def primes():
    found = []
    for n in count(2):
        if all((n % p != 0) for p in found):
            found.append(n)
            yield n

def neighbors(E):
    degs = primes()
    for l in degs:
        for f in E.isogenies_prime_degree(l):
            yield f

def is_dual(f, g):
    """A primitive probabilistic duality check, since computing the dual
    of an isogeny is complicated and expensive."""
    if f.degree() != g.degree():
        return False
    if f.domain() != g.codomain():
        return False
    if g.domain() != f.codomain():
        return False
    d = f.degree()
    E = f.domain()
    for _ in range(5):
        P = E.random_point()
        if P * d != g(f(P)):
            return False
    return True


def extract_cycle(bridge, path_mapping, start_E):
    """Assumes that path_mapping is a dictionary mapping j-invariants to isogenies
    whose codomain has this j-invariant. Furthermore, assumes that bridge is an
    isogeny such that both domain and codomain have entries in path_mapping, but
    the corresponding isogenies are different from bridge.

    In this case, the function traces the paths given by path_mapping and starting
    from bridge.domain() resp. bridge.codomain() backwards and combines them into
    a cycle in the isogeny graph (this assumes that tracing those paths backwards
    eventually stops at the same vertex). The returned isogenies are ordered such that
    the last isogeny in the list is the first edge of the cycle."""
    f = bridge
    E = f.domain()
    if path_mapping[f.codomain().j_invariant()] is not None:
        E2 = path_mapping[f.codomain().j_invariant()].codomain()
    else:
        E2 = start_E
    # first, we add isogenies E2 -> E -> ... -> start, using bridge as first isogeny
    f = E.isogeny(bridge.kernel_polynomial(), degree = bridge.degree(), codomain = E2)
    result = [f.dual()]
    current = E.j_invariant()
    while path_mapping[current] != None:
        result.append(path_mapping[current].dual())
        current = path_mapping[current].domain().j_invariant()
    result.reverse()
    # now add isogenies start -> ... -> E2, without using bridge 
    current = E2.j_invariant()
    while path_mapping[current] != None:
        result.append(path_mapping[current])
        current = path_mapping[current].domain().j_invariant()
    return result

def extract_path(E, path_mapping):
    """Assumes that path_mapping is a dictionary mapping j-invariants to isogenies
    whose codomain has this j-invariant.

    In this case, the function traces the path given by path mapping and starting
    at E backwards and gives a sequence of isogenies corresponding to this path."""
    f = path_mapping[E.j_invariant()]
    f = f.domain().isogeny(f.kernel_polynomial(), degree = f.degree(), codomain = E)
    result = [f]
    current = f.domain().j_invariant()
    while path_mapping[current] != None:
        result.append(path_mapping[current])
        current = path_mapping[current].domain().j_invariant()
    return result

def isogeny_graph_bfs(E, neighbor_iter = neighbors, find_cycle = False, target = None):
    if find_cycle == (target is not None):
        raise ValueError()
    if target is not None and not is_EllipticCurve(target):
        target = EllipticCurve(target.parent(), j = target)
    open = PriorityQueue()
    found = {}
    # add a random number as second parameter to break order ties in the first number;
    # this is necessary, otherwise PriorityQueue will try to compare EllipticCurves
    open.put((0, uniform(0, 1), (E, neighbor_iter(E))))
    found[E.j_invariant()] = None
    while True:
        prio, _, (current_E, edges) = open.get()
        f = next(edges)
        open.put((prio + 1, uniform(0, 1), (E, edges)))
        next_E = f.codomain()

        if next_E.j_invariant() in found:
            if found[current_E.j_invariant()] is not None and is_dual(f, found[current_E.j_invariant()]):
                continue
            if find_cycle:
                yield extract_cycle(f, found, E)
        elif target is not None and next_E.j_invariant() == target.j_invariant():
            found[next_E.j_invariant()] = f
            yield extract_path(target, found)
        else:
            found[next_E.j_invariant()] = f
            open.put((prio + 1, uniform(0, 1), (next_E, neighbor_iter(next_E))))

# p = 37
# F = GF(p**2)
# z = F.gen()
# j = 1 + 2*z
# E = EllipticCurve(F, j = j)
# E2 = EllipticCurve(F, j = j**p)
# path = isogeny_graph_bfs(E, target = E2)
# for f in reversed(path):
#     print(str(f.domain().j_invariant()) + " -> " + str(f.codomain().j_invariant()) + "; degree " + str(f.degree()))

