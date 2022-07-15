from modular_polynomial_sage import prime_modular_poly_mod_p

p = next_prime(20)
l = 2
e = int(log(p, l)) + 1
phi = prime_modular_poly_mod_p(l, p)
F = GF(p**2)
alpha = F.gen()
print(phi)
print(e)

P = PolynomialRing(F, [
    *['x' + str(i) for i in range(e)],
    *['z' + str(i) for i in range(e)], 
    *['y' + str(i) for i in range(1, e + 1)],
    *['w' + str(i) for i in range(1, e + 1)]
])
x = P.gens()[:e]
z = P.gens()[e:2*e]
y = [None, *P.gens()[2*e:3*e]]
w = [None, *P.gens()[3*e:]]

I = P.ideal([
    *[phi(x[i - 1] + z[i - 1] * alpha, x[i] + z[i] * alpha) for i in range(1, e)],
    phi(x[e - 1] + z[e - 1] * alpha, x[0] - z[0] * alpha),

    phi(x[0] + z[0] * alpha, y[1] + w[1] * alpha),
    *[phi(y[i] + w[i] * alpha, y[i + 1] + w[i + 1] * alpha) for i in range(1, e)],
    phi(y[e] + w[e] * alpha, x[0] - z[0] * alpha)
])
print(I.dimension())