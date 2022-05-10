
def division_polynomials(E, n):
    P = PolynomialRing(E.base_field(), 'x')
    x = P.gen(0)
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
        d = gcd(x1, gcd(y1_y, z1))
        (x1, y1_y, z1) = (P(x1/d), P(y1_y/d), P(z1/d))
        assert x1 in P and y1_y in P and z1 in P
        assert f(x1/z1) - y1_y**2 /f(x)/z1**2  == 0 
    return (x1, y1_y, z1)
