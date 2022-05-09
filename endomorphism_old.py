
class Endomorphism:

    def __init__(self, E, kernel_gens, insep_degree, sep_degree):
        assert E.is_ordinary()
        self.E = E
        self.kernel_gens = kernel_gens
        self.insep_degree = insep_degree
        self.sep_degree = sep_degree
        self.transforms = {}

    def div(self, n):
        P, Q = self.kernel_gens
        if P.order() % n != 0 or Q.order() % n != 0:
            raise ValueError
        return Endomorphism(self.E, (P * n, Q * n), self.insep_degree, self.sep_degree / n**2)

    def eval(self, point):
        if "phi1" not in self.__dict__:
            P, Q = self.kernel_gens
            F = P.curve().base_field()
            self._eval_phi1 = self.E.base_extend(F).isogeny(P)
            self._eval_E1 = self._eval_phi1.codomain()
            Q2 = self._eval_phi1(Q)
            self._eval_phi2 = self._eval_E1.isogeny(Q2)
            self._eval_E2 = self._eval_phi2.codomain()
            F = self.E.base_field()
            q = F.characteristic()**(log(self.insep_degree, F.characteristic()))
            self._eval_frobenius = lambda P: P.curve()(P[0]**q, P[1]**q, P[2]**q)
            self._eval_E3 = EllipticCurve(F, j = self._eval_E2.j_invariant()**q)
            self._eval_automorphism = sage.schemes.elliptic_curves.weierstrass_morphism.WeierstrassIsomorphism(E = self._eval_E3, F = self.E)
        return self._eval_automorphism(self._eval_frobenius(self._eval_phi2(self._eval_phi1(point))))

    def get_prime_power_transform(self, l, e):
        if l**e in self.transforms:
            return self.transforms[l**e]
        P, Q = torsion_group_gens(self.E, l, e)
        a, b = two_dim_discrete_log(self.eval(P), P, Q, l**e)
        c, d = two_dim_discrete_log(self.eval(Q), P, Q, l**e)
        R = Integers(l**e)
        self.transforms[l**e] = Matrix(R, [[a, c], [b, d]])
        return self.transforms[l**e]

    def deg(self):
        return self.sep_degree * self.insep_degree

    def add(self, rhs):
        # this is by Cauchy-Schwarz
        result_deg_bound = self.deg() + rhs.deg() + 2 * sqrt(self.deg() * rhs.deg())
        P0, Q0 = self.E(0, 1, 0), self.E(0, 1, 0)
        n0, m0 = 1, 1
        N = 1
        for (l, e) in prime_powers_leq(result_deg_bound):
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
        return Endomorphism(self.E, (P0, Q0), 0, N**2 / n0 / m0)