#![feature(fn_traits)]
#![feature(adt_const_params)]

extern crate feanor_la;

use feanor_la::prelude::*;
use feanor_la::wrapper::*;
use feanor_la::algebra::elliptic_curve::*;
use feanor_la::algebra::fq::*;
use feanor_la::algebra::poly::*;
use feanor_la::algebra::fractions::*;
use feanor_la::algebra::eea::*;

use feanor_la::la::vec::*;
use feanor_la::algebra::fq::zn_small::*;
use feanor_la::algebra::ring_ext::*;
use feanor_la::la::const_vector::*;

type P<K> = WrappingRing<PolyRing<K>>;

fn division_polynomials<'a, K: FiniteRing + CanonicalIsomorphismInfo<K>>(E: &EllipticCurve<WrappingRing<&'a K>>, n: usize) 
    -> (El<P<&'a K>>, El<P<&'a K>>, El<P<&'a K>>) 
{
    assert!(n > 0);
    let F = E.base_field();
    let P: P<&K> = PolyRing::adjoint(F.wrapped_ring().clone(), "X").bind_ring_by_value();
    let K = FieldOfFractions::new(P.wrapped_ring().clone()).bind_ring_by_value();
    let incl = embedding(F, &P);
    let x = P.wrapped_ring().bind_by_value(P.wrapped_ring().unknown());
    let A = incl(E.a4().clone());
    let B = incl(E.a6().clone());
    let f = x.pow(3) + &A * &x + &B;
    
    if n == 1 {
        return (x, f, P.one());
    }

    let mut x1 = x.pow(4) * 18 - &x * &f * 16 + x.pow(2) * &A * 12 + A.pow(2) * 2;
    let mut y1_y = -x.pow(6) * 27 - &A * x.pow(4) * 27 + x.pow(3) * &f * 28 - A.pow(2) * x.pow(2) * 9 + &A * &x * &f * 4 - A.pow(3) - &B * &f * 8;
    let mut z1 = &f * 8;

    for i in 2..n {
        let y1_squared = x1.pow(3) + &A * &x1 * z1.pow(2) + &B * z1.pow(3);
        println!("Found y1_squared");
        let x3 = (y1_squared - &y1_y * z1.pow(2) * 2 + z1.pow(3) * &f) - (&x1 - &z1 * &x).pow(2) * (&z1 * &x + &x1);
        println!("Found x3");
        let z3 = (&x1 - &z1 * &x).pow(2) * &z1;
        println!("Found z3");
        (x1, y1_y, z1) = (
            &x3 * (&x1 - &z1 * &x),
            ((&z1 * &f - &y1_y) * x3) + (&z3 * (y1_y * &x - &f * &x1)),
            z3 * (x1 - z1 * &x)
        );
        println!("Found (x1, y1_y, z1)");
        let d = gcd(&P, x1.clone(), gcd(&P, y1_y.clone(), z1.clone()));
        (x1, y1_y, z1) = (x1 / &d, y1_y / &d, z1 / &d);
        println!("Done step {}, degree is {}", i, x1.parent_ring().deg(x1.val()).unwrap_or(0));
    }
    return (x1, y1_y, z1)
}

type F37El = ZnElImpl<37, true>;
gen_const_vector!(ConstVector2F37; F37El; V0, V1);
type F1369MipoType = ConstVector2F37<{F37El::project(2)}, {F37El::project(33)}>;
const F1369_MIPO: Vector<F1369MipoType, F37El> = Vector::new(F1369MipoType {});
type F1369Type = SimpleRingExtension<StaticRing<F37El>, F1369MipoType, VectorArray<F37El, 2>>;
static F1369: F1369Type = F1369Type::new(F37El::RING, F1369_MIPO);

fn main() {
    let K = F1369.bind_ring();
    let a = F1369.bind(F1369.generator());
    let E = EllipticCurve::new(K, &a * 29 + 6, &a * 35 + 30);
    let (f, g, h): (El<P<&F1369Type>>, _, _) = division_polynomials(&E, 100);
    println!("{}", f.clone() / gcd(&f.ring(), f.clone(), h.clone()));
}
