read("modforms.gp");

E2(prec) = E_qexp(2,prec);

find_min_k(qexp,N,min_k=4,max_k=-1,profile=[5,25]) = {
    local(M,m,n,d,force,k,prec,bound);

    if(min_k < 4,error("***min_k has to be >= 4.***"));

    if(max_k >= 0,
        if(min_k > max_k,error("***min_k > max_k.***"),force=0);,
        force=1;
    );

    prec = poldegree(truncate(qexp),q)+1;
    qexp = Mod(qexp,N);
    k = min_k + (min_k%2); \\ Make sure k is even

    profile = concat(profile,[prec]);
    profile = vecsort(profile);

    while(k <= max_k || force,
        d = floor(k/12) + (k%12 != 2); \\ = dim M_k

        if(d>=prec,error("***Precision of the series is too low.***"));

        n = d;
        m = 1;
        while(m <= length(profile),
            bound = min(d + profile[m],prec);
            basis = victor_miller_basis(k,bound,N);
            while(n < bound && sum(i=0,d-1,polcoeff(qexp,i,q)*polcoeff(basis[i+1],n,q)) == polcoeff(qexp,n,q),
                n += 1;
            );
            if(n < bound, break(), if(n == prec, return(k)));
            m += 1;
        );
        if(n == prec, return(k));
        k += 2;
    );
    print("***No weight between ",min_k," and ",max_k," was found.***");
    return(-1);
}