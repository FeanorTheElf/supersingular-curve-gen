def class_group(O):
    assert type(O) == sage.rings.number_field.order.AbsoluteOrder_with_category
    assert O.discriminant() < 0
    K = O.number_field()
    assert K.degree() == 2

    h = O.class_number()

    a = 1
    b = K.maximal_order().ring_generator() * O.index_in(K.maximal_order())

    """
    The whole idea is to use Minkowski's theorem and observe that an integral ideal I
    contains an element a of norm at most

        (2/pi)^s [O_K : I] sqrt|d_K|

    For an ideal class c, find now I in c invertible with I^-1 integral. Then there is

        a in I with N(a) <= (2/pi)^s [O_K : I^-1] sqrt|d_K| = (2/pi)^s f N(I^-1) sqrt|d_K|

    Now note that I^-1 | (a) and so a/I^-1 = aI is integral.
    It has norm at most

        (2/pi)^s f sqrt|d_K| N(I^-1) N(I)

    The only problem now is that N(I^-1)N(I) does not have to be 1, as we are not necessarily
    in the maximal order.

    """

    class Ideal:

        def __init__(self, basis) -> None:
            self.basis = basis