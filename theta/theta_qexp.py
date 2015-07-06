"""This script can be called from the command line as follows:
1) With 3 arguments: > python theta_qexp.py M k N
This prints the q-expansions of all the theta series of parameter k attached
to the ideal classes of the fields of prime discriminant less than M in
absolute value. Only the first N coefficients are computed.

2) With 4 arguments: > python theta_qexp.py M k N file.txt
This command is similar to the previous one, except that the results are
printed in file.txt.
"""

import quadratic, primes
from math import sqrt

def coeff_form(f,k,N):
    """Return N coefficients of the theta series attached to the quad form.

    Given a quadratic form, return a list of length N, where the nth entry is
    the nth fourrier coefficient of the theta series corresponding to the
    form. See coeff_ideal for the details.

    One can also verifiy that using different representatives in the ideal
    class (i.e. choosing equivalent forms) gives proportional q-expansions.

    Examples:
    >>> c1 = coeff_form(quadratic.QuadraticForm(5, -6, 2), 2, 5)
    >>> c2 = coeff_form(quadratic.QuadraticForm(1, 0, 1), 2, 5)
    >>> c1; c2
    [0, -28+96*sqrt(-1), 112-384*sqrt(-1), 0, -448+1536*sqrt(-1)]
    [0, 4, -16, 0, 64]
    >>> x = quadratic.QuadraticRat(-4,-28, 96)
    >>> y = quadratic.QuadraticRat(-4, -448, 1536)
    >>> y/x
    16

    and note that 16 = 64/4.
    """

    return coeff_ideal(f.corresponding_ideal(), k, N)


def coeff_ideal(I,k,N):
    """Return N coefficients of the theta series attached to the ideal.

    Given an INTEGRAL ideal of the form I = [a,b + c*w_f], where w_f =f*w
    where w is the standard element of the integral basis of the corresponding
    quadratic integers, return a list of length N, where the nth entry is the
    sum of x^2k for all x in I such that N(x)/N(I) = n. See MScThesis for info
    on integral ideals.

    Note: The algorithm requires O(N) operations, where an operation is adding
    two powers of quadratic integers. This is when the ideal is fixed...

    Example:
    >>> I = quadratic.QuadraticIntId(D = -4, abc = (1,0,1))
    >>> coeff_ideal(I,0,10)
    [1, 4, 4, 0, 4, 8, 0, 0, 4, 4]
    >>> J = quadratic.QuadraticIntId(D = -3, abc = (1,0,1))
    >>> coeff_ideal(J,0,10)
    [1, 6, 0, 6, 6, 0, 0, 12, 0, 6]
    """

    # See ThetaSeriesqExp_comp2 for the computation of the bound
    a, b, c = I.get_abc()
    D = I.get_disc()
    f = I.corresponding_form()
    coeff = [quadratic.QuadraticOrder(D,0,0) for _ in range(N)]
    Bound = int(sqrt(2*N)/sqrt(f[0] + f[2] - sqrt((f[0] - f[2])**2 + f[1]**2)))
    for p in [(x,y) for x in range(-Bound,Bound + 1) for y in range(-Bound, Bound + 1)]:
        n = f(*p)
        if n < N:
            coeff[n] += quadratic.QuadraticOrder(D, a*p[0] - b*p[1], c*p[1])**(2*k)
    return coeff

def theta_expansions(D, k, N):
    """Compute the theta series q-expansions attached to the field of disc D.

    Given the discriminant of an order, one can compute the q-expansion of all
    the theta series attached to the ideal classes of this order. This
    function returns a string of the form

    'Theta series attached to the field of discriminant D:
        form1: [a0, a1, ..., aN-1]
        .
        .
        .
        formhD: [b0, b1, ... bN-1]
    '
    """

    res = ("Theta series (k = {k!s}) attached to the order of disc {D!s}\n"
           .format(k = k, D = D))
    for f in quadratic.QuadraticForm.Gaussian_forms(D):
        coeff = coeff_form(f, k, N)
        res += "    {:<15}: {coeff!s}\n".format(f.__str__(), coeff = coeff)
    return res



if __name__ == "__main__":
    import sys, doctest
    if len(sys.argv) >= 3 + 1: # M k N [name_of_file]
        M, k, N = list(map(int,sys.argv[1:4]))
        if len(sys.argv) == 5: # a file name is passed
            out = open(sys.argv[-1],'w')
        else:
            out = sys.stdout
        print("Printing results to",out.name)
        for D in [-p for p in primes.prime_list(M+1) if p % 4 == 3]:
            out.write(theta_expansions(D, k, N))
        out.close()
    doctest.testmod()
