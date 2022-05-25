
from numpy import poly


def extension_field(F, d):
    P = PolynomialRing(F, ['t'])
    t = P.gen(0)
    poly = t**d
    while not poly.is_prime():
        for i in range(min(d - 1, 10)):
            poly += F.random_element() * t**i
    return F.extension(poly, 'x')

def joint_field(F1, F2):
    F = F1.base_ring()
    assert F == F2.base_ring()
    if F1.degree() == 1:
        return (F2, F1.hom([F2(x) for x in F1.gens()], F2), F2.coerce_map_from(F2))
    f1 = F1.gen().minpoly()
    f2 = F2.gen().minpoly()
    assert gcd(f1.degree(), f2.degree()) == 1
    P = PolynomialRing(F, ['x', 'y'])
    x, y = P.gens()
    f1 = f1(x)
    f2 = f2(y)
    R = P.quo(P.ideal([f1, f2]))
    x, y = R(x), R(y)
    n = f1.degree() * f2.degree()
    A = Matrix(F, n + 1, n)
    while A.rank() < n - 1:
        potential_gen = F.random_element() * x + F.random_element() * y
        current = R(1)
        for power in range(n + 1):
            for i in range(f1.degree()):
                for j in range(f2.degree()):
                    A[power, i * f2.degree() + j] = P(current.lift()).coefficient({ x: i, y: j })
            current = current * potential_gen
    mipo_vec = kernel(A).basis()[0]
    x, = PolynomialRing(F, ['x']).gens()
    mipo = sum(x**i * mipo_vec[i] for i in range(n + 1))
    x_vec = vector(F, n)
    x_vec[f2.degree()] = 1
    x_im = A.solve_left(x_vec)
    y_vec = vector(F, n)
    y_vec[1] = 1
    y_im = A.solve_left(y_vec)
    result = F.extension(mipo, 'x')
    return (
        result,
        F1.hom([sum(x**i * x_im[i] for i in range(n + 1))], result),
        F2.hom([sum(x**i * y_im[i] for i in range(n + 1))], result)
    )

def joint_fields_multiple(F, degrees):
    embeddings = {}
    fields = {}
    current = extension_field(F, 1)
    for d in degrees:
        if d not in fields:
            new_field = extension_field(F, d)
            current, embedding1, embedding2 = joint_field(current, new_field)
            embeddings = { d: embedding1 * embeddings[d] for d in embeddings }
            embeddings[d] = embedding2
            fields[d] = new_field
    return (current, fields, embeddings)