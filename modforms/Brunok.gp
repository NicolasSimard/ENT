read("modforms.gp");

E2(prec) = E_qexp(2,prec);

find_min_k(qexp,N,min_k=4,max_k=-1,profile=[5,25]) = {
    local(X,m,n,d,force,k,prec,bound);

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
        \\if(k%500 == 0, print(k));

        d = floor(k/12) + (k%12 != 2); \\ = dim M_k

        if(d>=prec,print("***Precision of the series is too low.***"); return(-1));

        n = d;
        m = 1;
        while(m <= length(profile),
            bound = min(d + profile[m],prec);
            basis = victor_miller_basis(k,bound,N,0); \\ The basis is not reduced
            if(m == 1, X=matsolve(matrix(d,d,i,j,polcoeff(basis[j],i-1,q)),vector(d,i,polcoeff(qexp,i-1,q))~));
            while(n < bound && sum(i=0,d-1,X[i+1]*polcoeff(basis[i+1],n,q)) == polcoeff(qexp,n,q),
                n += 1;
            );
            if(n < bound, break(), if(n == prec, return(k)));
            m += 1;
        );
        if(n == prec, return(k));
        k += 2;
    );
    print("***No weight between ",min_k," and ",max_k," was found mod ",N,".***");
    return(-1);
}

