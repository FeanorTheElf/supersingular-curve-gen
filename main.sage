from itertools import count

class point(dict):
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

    def __init__(self, x, y):
        self['x'] = x
        self['y'] = y

p = 11

Fp = GF(p)
R = Fp['A', 'B', 'x', 'y']
A, B, x, y = R.gens()

def ec_add(phi, psi):
    if phi.x == psi.x and phi.y == psi.y:
        l = (3 * phi.x^2 + A) / (2 * phi.y)
    else:
        assert(phi.x != psi.x)
        l = (psi.y - phi.y) / (psi.x - phi.x)
    x = -phi.x - psi.x + l^2
    y = -phi.y + l * (phi.x - x)
    return point(x, y)

universal_point = point(x, y)

mul_p = universal_point
for _ in range(p - 1):
    mul_p = ec_add(mul_p, universal_point)
fp = point(x^p, y^p)
fp2 = point(x^(p^2), y^(p^2))
        
def elim_y(f):
    g = 0
    for k in count():
        coeff = f.coefficient({ y: k })
        if k % 2 == 0:
            g += coeff * (x^3 + A * x + B)^(k/2)
            f -= coeff * y^k
        else:
            assert(coeff == 0)
        if f == 0:
            return g

S = Fp['A', 'B', 'x']

def x_map(phi):
    """Returns an element in Frac(S) that describes the behavior of phi
    on the x-coordinate of points. This uniquely determines phi, as long
    as phi is an isogeny."""
    return S(elim_y(phi.x.numerator())) / S(elim_y(phi.x.denominator()))

U = S['T']
T, = U.gens()

mul_p_x = x_map(mul_p)
mipo = mul_p_x.numerator()(x = T) - mul_p_x * mul_p_x.denominator()(x = T)

print(x_map(ec_add(universal_point, universal_point)))
print(-2*x + (3*x^2 + A)^2/(x^3 + A*x + B)/4)

