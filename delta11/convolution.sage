def dir_convolution(L1,L2):
    """Compute the Dirichlet convolution of two Dirichlet L-series
    represented by their coefficients. Note that we need the nth
    coefficient to be at index n in the list. In particular, the
    lists must start with 0.

    Example:
    sage: dir_convolution([0]+[moebius(n) for n in range(1,8)],[0,1,1,1,1,1,1,1,1,1])
    
    [1, 0, 0, 0, 0, 0, 0]
    """
    if L2[0] != 0 or L2[0] != 0:
        return 0
    return [sum(L1[d]*L2[n//d] for d in divisors(n)) for n in range(1,min(len(L1),len(L2)))]
