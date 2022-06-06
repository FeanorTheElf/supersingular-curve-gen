#![feature(fn_traits)]
#![feature(adt_const_params)]
#![allow(non_snake_case)]

extern crate feanor_la;

pub mod cache;

use feanor_la::prelude::*;
use feanor_la::wrapper::*;
use feanor_la::elliptic_curve::*;
use feanor_la::ring_extension::simple_extension::*;
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

fn compute_factorization(E: EllipticCurve<WrappingRing<F10201Type>>) -> Vec<RingElWrapper<PolyRing<F10201Type>>> {
    let start = Instant::now();
    let ring = F10201.bind_ring_by_value();
    let z2 = ring.from(F10201.generator());
    let (mut f, _, h) = division_polynomials::division_polynomials(&E, 49);
    println!("Computed division polynomial in {} ms", start.elapsed().as_millis());
    let start = Instant::now();
    let d = gcd(&f.ring(), f.clone(), h.clone());
    f /= d;
    let ring = f.ring();
    let factorization = ring.factor(f);
    println!("Computed factorization in {} ms", start.elapsed().as_millis());
    for (g, _) in &factorization {
        println!("{}", g);
        println!("");
    }
    return factorization.into_iter().map(|(g, _)| g.into_val()).collect();
}

fn load_factorization() -> Vec<RingElWrapper<PolyRing<F10201Type>>> {
    cache::get_factorization()
}

fn main() {
    let ring = F10201.bind_ring_by_value();
    let z2 = ring.from(F10201.generator());
    let E = EllipticCurve::from_j_invariant(ring.clone(), z2 * 61 + 16);

    let factorization = load_factorization();
    let gen_poly = factorization.iter().filter(|p| p.deg().unwrap() == 21).next().unwrap();
    let splitting_field = SimpleRingExtension::<_, _, Vec<_>>::adjoin_element(F10201.clone(), |P| {
        P.embed(gen_poly.parent_ring(), gen_poly.val().clone())
    }, "z");
    let z = splitting_field.bind(splitting_field.generator());
    let i = embedding(ring.borrow_ring(), z.ring());
    let F = SimpleRingExtension::<_, _, VectorArray<_, 2>>::new(splitting_field.clone(), Vector::from_array([
        (-(z.pow(3) + i(E.a4().borrow_ring().clone()) * &z + i(E.a6().borrow_ring().clone()))).into_val(), 
        splitting_field.zero()
    ]), "w");
}
