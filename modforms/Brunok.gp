read("modforms.gp");


/*Given a prime p and a q-expansion f, return U_p(f). Note that this will
reduce the precision of the q-expansion.*/
U(p,f) = {
    local(prec);
    prec = poldegree(truncate(f),q)+1;
    return(Ser(vector(floor((prec-1)/p),n,polcoeff(f,(n-1)*p,q)),q));
}

V(p,f) = substpol(f,q,q^p);

E2(prec) = E_qexp(2,prec);

find_min_k(qexp,w=2,N,min_k=4,max_k=-1,profile=[5,25,-1]) = {
    local(phiN,X,m,n,d,force,k,prec,bound);

    if(min_k < 4,error("***min_k has to be >= 4.***"));

    if(max_k >= 0,
        if(min_k > max_k,error("***min_k > max_k.***"),force=0);,
        force=1;
    );

    prec = poldegree(truncate(qexp),q)+1;

    qexp = qexp*Mod(1,N);

    if(profile[length(profile)] < 0, profile[length(profile)] = prec);
    profile = vecsort(profile);

    phiN = eulerphi(N);
    k = ceil((min_k-w)/phiN)*phiN + w; \\ smallest number >= k_min and cong to
                                       \\ w mod phiN
    while(k <= max_k || force,
        print((k-w)/phiN);

        d = floor(k/12) + (k%12 != 2); \\ = dim M_k

        if(d>=prec,print("***No weight <",k,". Precision of the series is too low.***"); return(-1));

        n = d;
        m = 1;
        while(m <= length(profile),
            bound = min(d + profile[m],prec);
            basis = victor_miller_basis(k,bound,N,0); \\ The basis is not reduced
            if(m == 1, X=matsolve(matrix(d,d,i,j,polcoeff(basis[j],i-1,q)),vector(d,i,polcoeff(qexp,i-1,q))~));
            while(n < bound && sum(i=0,d-1,X[i+1]*polcoeff(basis[i+1],n,q)) == polcoeff(qexp,n,q),
                n += 1;
            );
            \\ At this point, n <= bound.
            if(n < bound,break(),
                if(n == prec || (m == length(profile) && n == bound), return([k,n-d]))
            );
            m += 1;
        );
        k += phiN;
    );
    print("***No weight between ",min_k," and ",max_k," was found mod ",N,".***");
    return(-1);
}

find_seq(qexp,w,p,M,known_seq=[],profile=[5,25,-1]) = {
    local(seq,n,k);
    if(M<=0,error("***M must be strictly positive. Recieved ",M,".***"));
    if(length(known_seq) == 0,
        seq = [find_min_k(qexp,w,p,4,max_k=-1,profile)];
        print("k_{",p,"^",1,"}=",seq[1]);,
        seq = known_seq;
    );
    for(n=length(seq)+1,M,
        k = find_min_k(qexp,w,p^n,seq[n-1][1],max_k=-1,profile);
        if(k[1] >= 4,
            seq = concat(seq,[k]);
            print("k_{",p,"^",n,"}=",seq[n]);,
            print("***The sequence could not be computed up to ",p,"^",M,". ***");
            print("   but only up to ",p,"^",n-1,". ***");
            return(seq);
        );
    );
    return(seq);
}
