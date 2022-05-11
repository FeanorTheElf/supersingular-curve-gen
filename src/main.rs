
extern crate feanor_la;

use feanor_la::prelude::*;
use feanor_la::wrapper::*;
use feanor_la::algebra::elliptic_curve::*;
use feanor_la::algebra::poly::*;
use feanor_la::algebra::fractions::*;

fn main() {
    let mut R = MultivariatePolyRing::new(i64::RING);
    let A = R.adjoint("A");
    let B = R.adjoint("B");
    let x1 = R.adjoint("x1");
    let x2 = R.adjoint("x2");
    let y1 = R.adjoint("y1");
    let y2 = R.adjoint("y2");
    let z1 = R.adjoint("z1");
    let z2 = R.adjoint("z2");
    let F1 = FieldOfFractions::new(R.clone());
    let e1 = embedding(&R, &F1);
    let F = F1.bind_ring_by_value();
    let e2 = F.wrapping_embedding();
    let e = compose::<_, _, MultivariatePolyRing<StaticRing<i64>>, FieldOfFractions<MultivariatePolyRing<StaticRing<i64>>>, WrappingRing<_>>(e2, e1);
    let (A, B, x1, x2, y1, y2, z1, z2) = (e(A), e(B), e(x1), e(x2), e(y1), e(y2), e(z1), e(z2));

    let E = EllipticCurve::new(&F, A, B);
    let P1 = EllipticCurvePoint::<&WrappingRing<_>>::Affine(x1/&z1, y1/z1);
    let P2 = EllipticCurvePoint::<&WrappingRing<_>>::Affine(x2/&z2, y2/z2);
    println!("{:?}", E.point_add(P1, P2));
<<<<<<< Updated upstream
    // ((y2^2z1^5z2^3 + -2 * y1y2z1^4z2^4 + y1^2z1^3z2^5 + -1 * x2^3z1^5z2^2 + x1x2^2z1^4z2^3 + x1^2x2z1^3z2^4 + -1 * x1^3z1^2z2^5) / (x2^2z1^5z2^3 + -2 * x1x2z1^4z2^4 + x1^2z1^3z2^5), 
    // (y2^3z1^9z2^4 + -3 * y1y2^2z1^8z2^5 + 3 * y1^2y2z1^7z2^6 + -1 * y1^3z1^6z2^7 + -1 * x2^3y2z1^9z2^3 + 2 * x2^3y1z1^8z2^4 + -3 * x1x2^2y1z1^7z2^5 + 3 * x1^2x2y2z1^7z2^5 + -2 * x1^3y2z1^6z2^6 + x1^3y1z1^5z2^7) / (-1 * x2^3z1^9z2^4 + 3 * x1x2^2z1^8z2^5 + -3 * x1^2x2z1^7z2^6 + x1^3z1^6z2^7))
=======
    // ((y2^2z1^5z2^3 + -2 * y1y2z1^4z2^4 + y1^2z1^3z2^5 + -1 * x2^3z1^5z2^2 + x1x2^2z1^4z2^3 + x1^2x2z1^3z2^4 + -1 * x1^3z1^2z2^5) / (x2^2z1^5z2^3 + -2 * x1x2z1^4z2^4 + x1^2z1^3z2^5), (y2^3z1^9z2^4 + -3 * y1y2^2z1^8z2^5 + 3 * y1^2y2z1^7z2^6 + -1 * y1^3z1^6z2^7 + -1 * x2^3y2z1^9z2^3 + 2 * x2^3y1z1^8z2^4 + -3 * x1x2^2y1z1^7z2^5 + 3 * x1^2x2y2z1^7z2^5 + -2 * x1^3y2z1^6z2^6 + x1^3y1z1^5z2^7) / (-1 * x2^3z1^9z2^4 + 3 * x1x2^2z1^8z2^5 + -3 * x1^2x2z1^7z2^6 + x1^3z1^6z2^7))
>>>>>>> Stashed changes
}
