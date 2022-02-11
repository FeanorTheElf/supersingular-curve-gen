
mod_polys = ClassicalModularPolynomialDatabase()

with open("output.log", "a") as file:

    d = 5
    p = next_prime(10000)
    file.write("d, p = %d, %d\n" % (d, p))
    R = PolynomialRing(GF(p), "x")
    x = R.gen()
    f = mod_polys[d](x, x^p)
    file.write(str(f) + "\n")

    R2 = PolynomialRing(GF(p^2), "x")
    f = R2(f)
    roots = f.roots()
    file.write(str(roots) + "\n")
    file.write("Using modular polynomial roots:\n")
    for (j, _) in roots:
        try:
            file.write(str(EllipticCurve(j = j).is_supersingular()) + "\n")
        except Exception as e:
            file.write(str(e) + "\n")

    file.write("Using field elements\n")
    for j in GF(p^2):
        file.write(str(EllipticCurve(j = j).is_supersingular()) + "\n")
