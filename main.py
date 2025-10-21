import gmpy2


RABIN_EXCEPTIONS = {
    5: (2, 1, 0, 0),
    10: (3, 1, 0, 0),
    13: (3, 2, 0, 0),
    34: (5, 3, 0, 0),
    37: (6, 1, 0, 0),
    58: (7, 3, 0, 0),
    61: (6, 5, 0, 0),
    85: (9, 2, 0, 0),
    130: (11, 3, 0, 0),
    214: (14, 3, 3, 0),
    226: (15, 1, 0, 0),
    370: (19, 3, 0, 0),
    526: (21, 9, 2, 0),
    706: (25, 9, 0, 0),
    730: (27, 1, 0, 0),
    829: (27, 10, 0, 0),
    1414: (33, 18, 1, 0),
    1549: (35, 18, 0, 0),
    1906: (41, 15, 0, 0),
    2986: (45, 31, 0, 0),
    7549: (85, 18, 0, 0),
    9634: (97, 15, 0, 0),
}


def euclid_filtered(a, b, bound):
    bound_sqrt = gmpy2.isqrt(bound)
    result = []
    while b != 0:
        r = a % b
        if r <= bound_sqrt:
            result.append(r)
        a, b = b, r
    return result


def two_squares(p):
    if p == 1:
        return [1, 0, 0, 0]
    elif p == 2:
        return [1, 1, 0, 0]

    rng = gmpy2.random_state(1)
    while True:
        q = gmpy2.mpz_random(rng, p)
        while gmpy2.powmod(q, (p - 1) // 2, p) != p - 1:
            q = gmpy2.mpz_random(rng, p)
        x = gmpy2.powmod(q, (p - 1) // 4, p)
        rems = euclid_filtered(p, x, p)
        if  len(rems) > 2:
            break

    return [rems[0], rems[1], 0, 0]

def three_squares(n):
    n = gmpy2.mpz(n)
    rng = gmpy2.random_state(1)
    if n in RABIN_EXCEPTIONS:
        return RABIN_EXCEPTIONS[n]
    elif n % 4 == 0:
        sub = three_squares(n // 4)
        return [sub[0] * 2, sub[1] * 2, sub[2] * 2, sub[3] * 2]
    elif n % 8 == 7:
        return [0, 0 , 0, 0]
    elif n % 8 == 3:
        
        while True:
            x = gmpy2.mpz_random(rng, gmpy2.isqrt(n) + 1)
            p = (n - x * x) // 2
            if not((n - x * x) % 2 != 0 or not (gmpy2.is_prime(p) or p == 1)):
                break
        two = two_squares(p)
        return [x, two[0] + two[1], abs(two[0] - two[1]), 0]
    elif gmpy2.is_square(n):
        return [gmpy2.isqrt(n), 0, 0, 0]
    else:
        while True:
            x = gmpy2.mpz_random(rng, gmpy2.isqrt(n) + 1)
            p = (n - x * x)
            if gmpy2.is_prime(p):
                break
        two = two_squares(p)
        return [x, two[0], two[1], 0]


def four_squares(n: int):
    n = gmpy2.mpz(n)

    if n == 1:
        return [1, 0, 0, 0]
    elif n == 2:
        return [1, 1, 0, 0]
    elif n == 3:
        return [1, 1, 1, 0]
    elif n == 0:
        return [0, 0, 0, 0]

    if n % 4 == 0:
        sub = four_squares(n // 4)
        return [*map(int, [sub[0] * 2, sub[1] * 2, sub[2] * 2, sub[3] * 2])]

    if gmpy2.is_square(n):
        return [*map(int, [gmpy2.isqrt(n), 0, 0, 0])]

    if n % 8 == 7:
        sub = three_squares(n - 4)
        return [*map(int, [sub[0], sub[1], sub[2], 2])]

    return [*map(int, three_squares(n))]
