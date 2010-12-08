def no_cancel_b_small(b, c, n, D, T):
    """
    Poly Risch Differential Equation - No cancelation: deg(b) small enough.

    Given a derivation D on k[t], n either an integer or +oo, and b, c
    in k[t] with deg(b) < deg(D) - 1 and either D == d/dt or
    deg(D) >= 2, either raise NonElementaryIntegral, in which case the
    equation Dq + b*q == c has no solution of degree at most n in k[t],
    or a solution q in k[t] of this equation with deg(q) <= n, or the
    tuple (h, b0, c0) such that h in k[t], b0, c0, in k, and for any
    solution q in k[t] of degree at most n of Dq + bq == c, y == q - h
    is a solution in k of Dy + b0*y == c0.
    """
    t = T[-1]
    d = D[-1]

    q = Poly(0, t)

    while not c.is_zero:
        if n == 0:
            m = 0
        else:
            m = c.degree(t) - d.degree(t) + 1

        if not 0 <= m <= n: # n < 0 or m < 0 or m > n
            raise NonElementaryIntegral

        if m > 0:
            p = Poly(c.as_poly(t).LC()/(m*d.as_poly(t).LC())*t**m, t)
        else:
            if b.degree(t) != c.degree(t):
                raise NonElementaryIntegral
            if b.degree(t) == 0:
                return (q, b.as_poly(T[-2]), c.as_poly(T[-2]))
            p = Poly(c.as_poly(t).LC()/b.as_poly(t).LC(), t)

        q = q + p
        n = m - 1
        c = c - derivation(p, D, T) - b*p

    return q
