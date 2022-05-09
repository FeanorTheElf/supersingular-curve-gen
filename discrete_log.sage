
def two_dim_discrete_log_baby_giant_step(P, Q1, Q2, n):
    """Assumes that Q1, Q2 are a basis of a free Z/nZ-module containing P.
    Computes a, b such that P = a Q1 + b Q2. This uses the baby-step-giant-step
    approach, so if n is not a prime, it will be faster to do it on each prime
    factor instead."""
    E = P.curve()
    assert P * n == E(0, 1, 0)
    assert Q1 * n == E(0, 1, 0)
    assert Q2 * n == E(0, 1, 0)
    steps = int(sqrt(n)) + 1
    R1 = Q1 * steps
    R2 = Q2 * steps
    giant_steps = { R1 * i + R2 * j: (i, j) for j in range(steps) for i in range(steps) }
    for k in range(steps):
        for l in range(steps):
            R = P + Q1 * k + Q2 * l
            if R in giant_steps:
                i, j = giant_steps[R]
                return (i * steps - k, j * steps - l)

def two_dim_discrete_log_prime_power(P, Q1, Q2, p, power):
    """Assumes that Q1, Q2 are a basis of a free Z/nZ-module containing P.
    Computes a, b such that P = a Q1 + b Q2. This assumes that n = p^power
    is a prime power, and will use power applications of the baby-step-giant-step
    method"""
    assert power >= 1
    if power == 1:
        return two_dim_discrete_log_baby_giant_step(P, Q1, Q2, p)

    # the basic idea is the exact sequence
    #   0 -> Z/pZ -> Z/p^(r+1)Z -> Z/p^rZ -> 0
    # where the two inner maps are multiplication by p^r and
    # mod p^r, respectively
    right_P, right_Q1, right_Q2 = P * p, Q1 * p, Q2 * p
    right_a, right_b = two_dim_discrete_log_prime_power(right_P, right_Q1, right_Q2, p, power - 1)
    left_P = P - right_a * Q1 - right_b * Q2
    # left_P is in the kernel of Z/p^(r+1)Z -> Z/p^rZ, thus in the image of Z/pZ -> Z/p^(r+1)Z
    left_Q1 = Q1 * p**(power - 1)
    left_Q2 = Q2 * p**(power - 1)
    left_a, left_b = two_dim_discrete_log_baby_giant_step(left_P, left_Q1, left_Q2, p)
    return (left_a * p**(power - 1) + right_a, left_b * p**(power - 1) + right_b)

def two_dim_discrete_log(P, Q1, Q2, n):
    """Assumes that Q1, Q2 are a basis of a free Z/nZ-module containing P.
    Computes a, b such that P = a Q1 + b Q2"""
    a, b = (0, 0)
    m = 1
    for (p, power) in factor(n):
        k = p**power
        a1, b1 = two_dim_discrete_log_prime_power(P * Integer(n/k), Q1 * Integer(n/k), Q2 * Integer(n/k), p, power)
        a = crt(a, a1, m, k)
        b = crt(b, b1, m, k)
        m *= k
    assert m == n
    return (a, b)