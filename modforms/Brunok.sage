def GG(k,prec):
    """ Returns the coefficients of the q-expansion of the Eisenstein series
    of weight k up to a given precision, where the series is normalized to
    be a newform (i.e. the coefficient of q is 1).

    INPUT:
    - k     : the weight
    - prec  : the number of coefficients

    OUTPUT:
    - L     : list of the coefficients (all integral except the first)
    """
    c0 = -1/24 if k==2 else -bernoulli(k)/(2*k)
    return [c0] + [sigma(n,k-1) for n in range(1,prec)]

def Eisen_coeff(k,prec):
    """ Returns the coefficients of the q-expansion of the Eisentein series
    of weight k up to a given precision, where the series is now normalized
    to have constant term = 1 (then all the coefficients are inetgral).

    INPUT:
    - k     : the weight
    - prec  : the number of coefficients

    OUTPUT:
    - L     : list of the coefficients (all integral)
    """
    c0 = GG(k,prec)[0]
    return [c/c0 for c in GG(k,prec)]


def find_min_k(coeff,N,k=4,max_k=200):
    """This function returns a modular form congruent mod N to
    a given modular form (or q-expansion, more generally).

    INPUT:
    - coeff : coefficients of the modular form (list of rationnal nbrs)
    - N     : the modulus
    - k (optional): starting weight (default =4)
    - max_k (optional): maximum weight that we want (default = 200)

    OUTPUT:
    - (k,f) : a genuine modular form f of weight k congruent to coeff mod N

    TEST:
    sage: delta=CuspForms(1,12).basis()[0].qexp(10).list()
    sage: find_min_k(delta,100)
    (12,
     q - 24*q^2 + 252*q^3 - 1472*q^4 + 4830*q^5 - 6048*q^6 - 16744*q^7 + 84480*q^8 - 113643*q^9 + O(q^10))
    sage: f = find_min_k([1,2,3,4,5],3)[1].list();
    sage: are_equal_mod(f,[1,2,3,4,5],3)
    True
    """

    while k <= max_k:
        if k%2 != 0: continue
        # Precision of the given q-expansion
        prec = len(coeff)
        # Find the basis up to O(q^prec). This gives prec coefficients
        basis = victor_miller_basis(k,prec)
        # Dimension of the space of modular forms of level 1, weight k
        d = len(basis)
        # If the dimension is larger than the number of coefficients, we are
        # certain to find an answer (due to the nature of the Victor Miller
        # basis), but this might not be the right answer!
        if d >= prec:
            print "***Warning: precision of the series is low; Answer might be meaningless.***"
        # Extract the coefficients of the basis elements
        basis_coeffs = map(lambda f: f.list(),basis)
        # Check if the coefficients match mod N up to the precision
        for n in range(d,prec):
            # If they don't, we stop the loop
            if sum(coeff[i]*basis_coeffs[i][n] for i in range(d))%N != coeff[n]%N:
                break
        else:
            # If they match up to the precision, we return k and the modular form
            # which reduces to the given form mod N
            return (k, sum(coeff[i]*basis[i] for i in range(d)))
        # At this point in the while loop, we know that k is even
        k = k+2
    # At this point, k reached the mak_k bound without finding anything...
    print "No weight smaller than " + str(max_k) + " was found..."

def are_equal_mod(L1,L2,N):
    """ Returns true if the two lists are congruents mod N.
    """
    for t in zip(L1,L2):
        if t[0]%N != t[1]%N:
            return False
    else:
        return True
