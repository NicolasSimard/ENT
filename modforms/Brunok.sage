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

def E2(prec):
    return Eisen_coeff(2,prec)

def find_min_k(coeff,N,min_k=4,max_k=None,output_lift=False):
    """This function returns a modular form congruent mod N to
    a given modular form (or q-expansion, more generally).

    INPUT:
    - coeff        : coefficients of the modular form (list of rationnal nbrs)
    - N            : the modulus
    - min_k (opt.) : starting weight (default =4)
    - max_k (opt.) : maximum weight that we want (default = None)
    - output_lift (opt.): return the form that reduces to the input form mod N.

    OUTPUT:
    - k : a weight k s.t. there is a modular form of weight k congruent to coeff mod N
          If output_lift = True, the output is (k,f), where f is the modular form in question.

    EXAMPLES:
    1) To find a weight at all cost. Be careful if you are not sure this weight exists:
    find_min_k(coeff,5^3)
    
    2) To find a weight between 4 and 200:
    find_min_k(coeff,5^3,max_k=200)

    3) To find a weight between 100 and 200:
    find_min_k(coeff,5^3,min_k=100,max_k=200) 
    
    4) To find a weight greater than 200 at all cost:
    find_min_k(coeff,5^3,min_k=200)

    TEST:
    sage: delta=CuspForms(1,12).basis()[0].qexp(10).list()
    sage: find_min_k(delta,100)
    (12,
     q - 24*q^2 + 252*q^3 - 1472*q^4 + 4830*q^5 - 6048*q^6 - 16744*q^7 + 84480*q^8 - 113643*q^9 + O(q^10))
    sage: f = find_min_k([1,2,3,4,5],3,output_lift=True)[1].list();
    sage: are_equal_mod(f,[1,2,3,4,5],3)
    True
    """

    if min_k < 4:
        raise ValueError("***min_k has to be >= 4.***")

    if max_k is not None: # k_max has been set
        if min_k > max_k: # It has been set to an invalid value
            raise ValueError("***min_k > max_k.***")
        force = False
    else: # max_k has not been set, so we force computations
        force = True

    coeffs_modN = [Mod(c,N) for c in coeff]

    k = min_k
    while force or k <= max_k:
        if k%2 != 0: 
            k = k+1
            continue

        # Precision of the given q-expansion
        prec = len(coeff)
        # Dimension of the space of modular forms of level 1, weight k
        d = k//12 + (1 if k%12 != 2 else 0)
        # If the dimension is larger than the number of coefficients, we are
        # certain to find an answer (due to the nature of the Victor Miller
        # basis), but this might not be the right answer!
        if d >= prec:
            print "***Precision of the series is too low to find k.***"
            break
        coeff_part = [(d,min(d+5,prec)),(min(d+5,prec),min(d+25,prec)),(min(d+25,prec),prec)]
        # coeff_part = [(d,d+5),(d+5,prec)] # Seems slower (empirically)
        n = d - 1
        for part in coeff_part:
            basis = victor_miller_basis(k,part[1])
            basis_modN = [f*Mod(1,N) for f in basis]
            while n < part[1] and sum(coeffs_modN[i]*basis_modN[i][n] for i in xrange(d)) == coeffs_modN[n]:
                n = n + 1
            if n < part[1]:
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
    print "***No weight between " + str(min_k) + " and " + str(max_k) + " was found mod " + str(N) + ".***"
    return None

def find_seq(coeff,p,M,max_k=None,known_seq=[],verb=False):
    """ Given the coefficients of the q-expansion of a p-adic modular form, returns
    the sequence of minimal weights for modulus p, p^2, ..., p^M.

    INPUT:
    - coeff       : Coefficients of the modular form.
    - p           : prime number.
    - M           : M such that modulus goes up to p^M.
    - max_k (opt.): the maximal weight that the algorithm will reach.
    - known_seq (opt.): Previously computed sequence. Saves time!
    - verb (opt.) : print the computed values of the sequence along the way. (default = false)

    OUTPUT:
    - seq         : Sequence of integers [kn] such that kn is the minimal weight mod p^(n+1)

    sage: find_seq(Eisen_coeff(2,200),3,5)
    [4, 8, 20, 56, 164]
    sage: find_seq(Eisen_coeff(2,200),5,5)
    [6, 42, 202, 1002, 2014]
    """

    if len(known_seq) == 0:
        seq = [find_min_k(coeff,p,max_k=max_k)]
        if verb:
            print "k_{" + str(p) + "^" + str(1) + "}: " + str(seq[-1])
    else:
        seq = known_seq
    for n in range(len(seq)+1,M+1):
        k=find_min_k(coeff,p^n,min_k=seq[-1],max_k=max_k)
        if k is not None:
            seq.append(k)
            if verb:
                 print "k_{" + str(p) + "^" + str(n) + "}: " + str(seq[-1])
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
