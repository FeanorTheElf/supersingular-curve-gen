#![feature(fn_traits)]
#![feature(adt_const_params)]
#![feature(unboxed_closures)]
#![allow(non_snake_case)]

extern crate feanor_la;

use feanor_la::prelude::*;
use feanor_la::wrapper::*;
use feanor_la::fq::FiniteRing;
use feanor_la::fq::zn_big::*;
use feanor_la::elliptic_curve::*;
use feanor_la::ring_extension::simple_extension::*;
use feanor_la::ring_extension::extension_wrapper::*;
use feanor_la::discrete_log::*;
use feanor_la::finite_field_sqrt::sqrt;
use feanor_la::poly::*;
use feanor_la::eea::*;

use std::time::Instant;

pub mod f10201 {
    use feanor_la::fq::fq_small::define_fq::*;

    type F101El = ZnElImpl<101, true>;
    gen_const_vector!(ConstVector2F101; F101El; V0, V1);
    type F10201MipoType = ConstVector2F101<{F101El::project(-2)}, {F101El::project(-97)}>;
    const F10201_MIPO: Vector<F10201MipoType, F101El> = Vector::new(F10201MipoType::INSTANCE);
    pub type F10201Type = SimpleRingExtension<StaticRing<F101El>, F10201MipoType, VectorArray<F101El, 2>>;
    pub static F10201: F10201Type = F10201Type::new(F101El::RING, F10201_MIPO, "α");
    pub static F101: StaticRing<F101El> = F101El::RING;
}

use f10201::{F10201Type, F10201};

type BaseField = WrappingRing<&'static F10201Type>;
type PolyRingBaseField = WrappingRing<PolyRing<&'static F10201Type>>;

fn compute_division_poly_factorization(E: &EllipticCurve<BaseField>, n: usize) -> Vec<El<PolyRingBaseField>> {
    let start = Instant::now();
    let (f, _, mut h) = division_polynomials::division_polynomials(&E, n);
    println!("Computed division polynomial in {} ms", start.elapsed().as_millis());
    let start = Instant::now();
    let d = gcd(&f.ring(), f, h.clone());
    h /= d;
    let ring = h.ring();
    let factorization = ring.factor(h);
    println!("Computed factorization in {} ms", start.elapsed().as_millis());
    return factorization.into_iter().map(|(g, _)| g.into_val()).collect();
}

fn compute_frobenius_transform(E: &EllipticCurve<BaseField>, p: usize, e: u32) -> Matrix<MatrixOwned<El<Zn>>, El<Zn>> { 
    let factorization = compute_division_poly_factorization(&E, p.pow(e));
    let deg = factorization.iter()
        .filter_map(|f| f.deg())
        .filter(|d| *d > 0)
        .fold(1, |a, b| lcm(&i64::RING, a, b as i64))
        .abs() as usize;
    let mut gen_poly = factorization.iter().filter(|p| p.deg().unwrap() == deg).next().unwrap().clone();
    println!("Generator polynomial: {}", gen_poly);
    let scale_factor = gen_poly.lc().unwrap().inv();
    gen_poly = gen_poly.scaled(&scale_factor);
    let ext_field = ExtensionWrapper::from(F10201.clone()).extend(
        |ring| SimpleRingExtension::adjoin_element(ring, |_| gen_poly.into_val(), "t")
    );
    let point_field = ext_field.bind_ring();
    match compute_torsion_basis_in_field(&E, &factorization, &point_field, p as i64, e) {
        Ok((P1, P2)) => return compute_frobenius_transform_of_basis(&E, P1, P2, &point_field, p as i64, e),
        Err(irred_poly) => {
            println!("Cannot take root of {}, going to field extension", irred_poly);
            assert_eq!(2, irred_poly.deg().unwrap());
            let poly = irred_poly.into_val();
            let ext_field = ext_field.extend(
                |ring| SimpleRingExtension::adjoin_element(ring, |_| poly, "w")
            );
            let point_field = ext_field.bind_ring();
            let (P1, P2) = compute_torsion_basis_in_field(&E, &factorization, &point_field, p as i64, e).unwrap();
            return compute_frobenius_transform_of_basis(&E, P1, P2, &point_field, p as i64, e);
        }
    };
}

fn compute_torsion_basis_in_field<'a, L>(
    E: &EllipticCurve<BaseField>, 
    factorization: &Vec<El<PolyRingBaseField>>,
    ext_field: &WrappingRing<&'a L>, 
    p: i64, 
    e: u32,
) -> Result<(EllipticCurvePoint<WrappingRing<&'a L>>, EllipticCurvePoint<WrappingRing<&'a L>>), El<WrappingRing<PolyRing<&'a L>>>>
    where L: FiniteRing + CanonicalEmbeddingInfo<F10201Type> + HashableElRing + DivisibilityInfoRing
{
    let poly_ring = PolyRing::adjoint(*ext_field.wrapped_ring(), "X").bind_ring_by_value();
    let i = compose::<_, _, BaseField, WrappingRing<&L>, WrappingRing<PolyRing<&L>>>(
        embedding(ext_field, &poly_ring),
        embedding(E.base_field(), ext_field)
    );
    let lift = poly_ring.lift_hom::<PolyRing<&F10201Type>, _>(
        poly_ring.wrapped_ring().lift_hom::<_, &F10201Type>(
            embedding(&F10201, *ext_field.wrapped_ring())
        )
    );
    let X = poly_ring.from(poly_ring.wrapped_ring().unknown());
    let E_poly: El<WrappingRing<PolyRing<&L>>> = X.pow(3) + i(E.a4().clone()) * &X + i(E.a6().clone());

    let mut points = factorization.iter()
        .flat_map(|f| lift(f.clone()).roots())
        .map(|(r, _)| r)
        .map(|x| sqrt(E_poly(x.clone()), &ext_field).map(|y| EllipticCurvePoint::Affine(x.clone(), y)).ok_or(X.pow(2) - poly_ring.embed(ext_field, E_poly(x.clone()))))
        .inspect(|P| if let Ok(P) = P { assert!(E.is_on_curve(&P, &ext_field)) });
    
    let n = BigInt::from(p.pow(e - 1));
    let mut P1 = points.next().unwrap()?;
    while E.mul_point(&P1, &n, &ext_field) == EllipticCurvePoint::Infinity {
        P1 = points.next().unwrap()?;
    }
    let mut P2 = points.next().unwrap()?;
    while E.mul_point(&P2, &n, &ext_field) == EllipticCurvePoint::Infinity || discrete_log(
        E.mul_point(&P2, &n, &ext_field), 
        &P1, 
        &i64::RING.bind(p.pow(e) as i64), 
        |a, b| E.point_add(a, b, &ext_field), 
        EllipticCurvePoint::Infinity
    ).is_some() {
        P2 = points.next().unwrap()?;
    }
    return Ok((P1, P2));
}

fn frobenius<K>(P: &EllipticCurvePoint<WrappingRing<K>>, q: &BigInt) -> EllipticCurvePoint<WrappingRing<K>>
    where K: FiniteRing
{
    match P {
        EllipticCurvePoint::Infinity => EllipticCurvePoint::Infinity,
        EllipticCurvePoint::Affine(x, y) => EllipticCurvePoint::Affine(x.pow_big(q), y.pow_big(q))
    }
}

fn compute_frobenius_transform_of_basis<'a, L>(
    E: &EllipticCurve<BaseField>, 
    P1: EllipticCurvePoint<WrappingRing<&'a L>>, 
    P2: EllipticCurvePoint<WrappingRing<&'a L>>, 
    field: &WrappingRing<&'a L>, 
    p: i64, 
    e: u32
) -> Matrix<MatrixOwned<El<Zn>>, El<Zn>>
    where L: FiniteRing + CanonicalEmbeddingInfo<F10201Type> + HashableElRing
{
    assert_eq!(E.mul_point(&P1, &BigInt::from(p.pow(e)), &field), EllipticCurvePoint::Infinity);
    assert_eq!(E.mul_point(&P2, &BigInt::from(p.pow(e)), &field), EllipticCurvePoint::Infinity);
    assert_ne!(E.mul_point(&P1, &BigInt::from(p.pow(e - 1)), &field), EllipticCurvePoint::Infinity);
    assert_ne!(E.mul_point(&P2, &BigInt::from(p.pow(e - 1)), &field), EllipticCurvePoint::Infinity);

    let q = F10201.size();
    let Q1 = frobenius(&P1, &q);
    let Q2 = frobenius(&P2, &q);

    let (a, b) = discrete_log_2d(
        Q1, 
        (&P1, &P2), 
        (&i64::RING.bind_by_value(p.pow(e)), &i64::RING.bind_by_value(p.pow(e))), 
        |a, b| E.point_add(a, b, &field), 
        EllipticCurvePoint::Infinity
    ).unwrap();

    let (c, d) = discrete_log_2d(
        Q2, 
        (&P1, &P2), 
        (&i64::RING.bind_by_value(p.pow(e)), &i64::RING.bind_by_value(p.pow(e))), 
        |a, b| E.point_add(a, b, &field), 
        EllipticCurvePoint::Infinity
    ).unwrap();

    let ring = Zn::new(BigInt::from(p.pow(e)));
    let project = |x: RingElWrapper<StaticRing<i64>>| ring.from_z(*x.val());

    return Matrix::from_array([[project(a), project(c)], [project(b), project(d)]]);
}

fn load_torsion_basis() -> (
    EllipticCurvePoint<SimpleRingExtension<F10201Type, Vec<El<F10201Type>>, Vec<El<F10201Type>>>>, 
    EllipticCurvePoint<SimpleRingExtension<F10201Type, Vec<El<F10201Type>>, Vec<El<F10201Type>>>>,
    SimpleRingExtension<F10201Type, Vec<El<F10201Type>>, Vec<El<F10201Type>>>
) {
    let poly_ring: PolyRingBaseField = PolyRing::adjoint(&F10201, "x").bind_ring_by_value();
    let x = poly_ring.from(poly_ring.wrapped_ring().unknown());
    let a = poly_ring.from(poly_ring.wrapped_ring().from(F10201.generator()));
    let poly = (&a * 96 + 21) + (&a * 78 + 18) * &x + (&a * 50 + 83) * x.pow(2) + (&a * 42 + 69) * x.pow(3) + (&a * 15 + 88) * x.pow(4) + (&a * 48 + 20) * x.pow(5) + (&a * 54 + 82) * x.pow(6) + (&a * 95 + 4) * x.pow(7) + (&a * 16 + 52) * x.pow(8) + (&a * 18 + 15) * x.pow(9) + (&a * 13 + 76) * x.pow(10) + (&a * 93 + 84) * x.pow(11) + (&a * 76 + 24) * x.pow(12) + (&a * 51 + 87) * x.pow(13) + (&a * 17 + 27) * x.pow(14) + (&a * 3 + 68) * x.pow(15) + (&a * 33 + 75) * x.pow(16) + (&a * 24 + 86) * x.pow(17) + (&a * 38 + 85) * x.pow(18) + (&a * 26 + 91) * x.pow(19) + (&a * 92 + 61) * x.pow(20) + x.pow(21);
    
    let field: SimpleRingExtension<F10201Type, Vec<_>, Vec<_>> = SimpleRingExtension::adjoin_element(
        F10201.clone(), 
        |ring| ring.lift_hom::<_, &F10201Type>(embedding(&F10201, &F10201))(poly.val().clone()), 
        "x"
    );
    let extension_field = field.bind_ring();
    let x = extension_field.from(extension_field.wrapped_ring().generator());
    let a = extension_field.from(extension_field.wrapped_ring().from(F10201.generator()));
    return (
        EllipticCurvePoint::Affine(
            x.clone().into_val(), 
            ((&a * 38 + 19) + (&a * 32 + 98) * &x + (&a * 11 + 20) * x.pow(2) + (&a * 36 + 62) * x.pow(3) + (&a * 86 + 72) * x.pow(4) + (&a * 94 + 67) * x.pow(5) + (&a * 49 + 4) * x.pow(6) + (&a * 79 + 64) * x.pow(7) + (&a * 10 + 47) * x.pow(8) + (&a * 57 + 91) * x.pow(9) + (&a * 80 + 10) * x.pow(10) + (&a * 40 + 100) * x.pow(11) + (&a * 76 + 65) * x.pow(12) + x.pow(13) * 38 + (&a * 67 + 15) * x.pow(14) + (&a * 75 + 62) * x.pow(15) + (&a * 70 + 38) * x.pow(16) + (&a * 97 + 75) * x.pow(17) + (&a * 26 + 90) * x.pow(18) + (&a * 91 + 12) * x.pow(19) + (&a * 20 + 32) * x.pow(20)).into_val()
        ),
        EllipticCurvePoint::Affine(
            ((&a * 57 + 8) + (&a * 28 + 97) * &x + (&a * 12 + 79) * x.pow(2) + (&a * 24 + 35) * x.pow(3) + (&a * 52 + 9) * x.pow(4) + (&a * 100 + 97) * x.pow(5) + (&a * 60 + 48) * x.pow(6) + (&a * 71 + 57) * x.pow(7) + (&a * 84 + 85) * x.pow(8) + (&a * 6 + 30) * x.pow(9) + (&a * 9 + 55) * x.pow(10) + (&a * 84 + 96) * x.pow(11) + (&a * 62 + 24) * x.pow(12) + (&a * 83 + 100) * x.pow(13) + (&a * 94 + 100) * x.pow(14) + (&a * 50 + 20) * x.pow(15) + (&a * 84 + 54) * x.pow(16) + (&a * 11 + 97) * x.pow(17) + (&a) * 26 * x.pow(18) + (&a * 50 + 98) * x.pow(19) + (&a * 4 + 100) * x.pow(20)).into_val(), 
            ((&a * 60 + 3) + (&a * 23 + 42) * &x + (&a * 47 + 77) * x.pow(2) + (&a * 70 + 66) * x.pow(3) + (&a) * 30 * x.pow(4) + (&a * 84 + 20) * x.pow(5) + (&a * 14 + 90) * x.pow(6) + (&a * 86 + 53) * x.pow(7) + (&a * 46 + 43) * x.pow(8) + (&a * 65 + 81) * x.pow(9) + (&a * 76 + 21) * x.pow(10) + (&a * 87 + 83) * x.pow(11) + (&a * 48 + 44) * x.pow(12) + (&a * 33 + 1) * x.pow(13) + (&a * 60 + 20) * x.pow(14) + (&a * 37 + 6) * x.pow(15) + (&a * 55 + 75) * x.pow(16) + (&a * 16 + 31) * x.pow(17) + (&a * 51 + 49) * x.pow(18) + (&a * 24 + 50) * x.pow(19) + (&a * 14 + 30) * x.pow(20)).into_val()
        ),
        (**extension_field.wrapped_ring()).clone()
    );
}

fn is_square(x: i64) -> bool {
    let root = ((x as f64).sqrt() as i64);
    return root * root == x || (root + 1) * (root + 1) == x;
}

fn is_nontrivially_solvable(D: i64, m: i64) -> bool {
    (0..=((m as f64).sqrt() as i64)).map(|x| m - x * x).filter(|Dy| *Dy != 0).any(|Dy| Dy % D == 0 && is_square(Dy / D))
}

fn has_small_l_vulcano(D: i64, n: i64, l: i64) -> bool {
    let max_crater_size = (n as f64).ln() as i64 + 1;
    for e in 3..max_crater_size {
        if is_nontrivially_solvable(D, l.pow(e as u32)) {
            return true;
        }
    }
    return false;
}

fn main() {
    for n in (30..100000).step_by(100) {
        let mut n = n;
        while !i64::RING.is_prime(&(n as i64)) {
            n += 1;
        }
            let mut l = (n as f64).ln() as i64;
            while !i64::RING.is_prime(&l) {
                l += 1;
            }
            let c = (1..=n.pow(2)).filter(|D| has_small_l_vulcano(*D, n, l)).count();
            println!("{}, {}, {}", n, (c as f64) / (n as f64).sqrt(), c);
    }
}
