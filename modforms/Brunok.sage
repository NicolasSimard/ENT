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
    to have constant term = 1.

    INPUT:
    - k     : the weight
    - prec  : the number of coefficients

    OUTPUT:
    - L     : list of the coefficients (all integral)
    """
    c0 = GG(k,prec)[0]
    return [c/c0 for c in GG(k,prec)]


def find_min_k(coeff,N,k=4,max_k=200,force = False,output_lift=False):
    """This function returns a modular form congruent mod N to
    a given modular form (or q-expansion, more generally).

    INPUT:
    - coeff        : coefficients of the modular form (list of rationnal nbrs)
    - N            : the modulus
    - k (opt.)     : starting weight (default =4)
    - max_k (opt.) : maximum weight that we want (default = 200)
    - force (opt.) : If set to True, will force the algorithm to find k.
    - output_lift (opt.): return the form that reduces to the input form mod N.

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

    while k <= max_k or force:
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
        basis_coeffs = [f*Mod(1,N) for f in basis]
        coeffs_modN = [Mod(c,N) for c in coeff]
        # Check if the coefficients match mod N up to the precision
        for n in range(d,prec):
            # If they don't, we stop the loop
            if sum(coeffs_modN[i]*basis_coeffs[i][n] for i in xrange(d)) != coeffs_modN[n]:
                break
        else:
            # If they match up to the precision, we return k and the modular form
            # which reduces to the given form mod N
            if output_lift:
                return (k, sum(coeff[i]*basis[i] for i in range(d)))
            else:
                return k
        # At this point in the while loop, we know that k is even
        k = k+2
    # At this point, k reached the mak_k bound without finding anything...
    print "***No weight smaller than " + str(max_k) + " was found mod" + str(N) + ".***"
    return None

def find_seq(coeff,p,M,max_k=200,known_seq=[],force=False):
    """ Given the coefficients of the q-expansion of a p-adic modular form, returns
    the sequence of minimal weights for modulus p, p^2, ..., p^M.

    INPUT:
    - coeff       : Coefficients of the modular form
    - p           : prime number
    - M           : M such that modulus goes up to p^M
    - max_k (opt.): the maximal weight that the algorithm will reach
    - known_seq (opt.): Previously computed sequence. Saves time!
    - force (opt.): will force find_min_k (see above)

    OUTPUT:
    - seq         : Sequence of integers [kn] such that kn is the minimal weight mod p^(n+1)

    sage: find_seq(Eisen_coeff(2,200),3,7,force=True)
    [4, 8, 20, 56, 164, 488, 974]
    sage: find_seq(Eisen_coeff(2,200),5,7,force=True)
    [6, 42, 202, 1002, 2014, 2286, 2350]
    """
    if len(known_seq) == 0:
        seq = [find_min_k(coeff,p,max_k=max_k,force=force)[0]]
    else:
        seq = known_seq
    for n in range(len(seq)+1,M+1):
        min_k=find_min_k(coeff,p^n,k=seq[-1],max_k=max_k,force=force)
        if min_k is not None:
            seq.append(min_k[0])
        else:
            print "***Error: the sequence could not be computed up to " + str(p) + "^" + str(M) + ","
            print "but only up to " + str(p) + "^" + str(n-1) + ". Consider increasing max_k.***"
            return seq
    return seq

def are_equal_mod(L1,L2,N):
    """ Returns true if the two lists are congruents mod N.
    """
    for t in zip(L1,L2):
        if t[0]%N != t[1]%N:
            return False
    else:
        return True
