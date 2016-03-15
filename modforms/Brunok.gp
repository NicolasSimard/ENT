\r vm_basis.gp
\r operators.gp

E2star(p) = mfadd(E(2),mfmul(-p,V(p)(E(2))));
pjVjE2star(p,j) = mfmul(p^j,V(p,j)(E2star(p)));

find_min_k(f,w,N,match='auto,k_range=[4,[1]]) =
{
    my(phiN,X,m,n,d,pr,k,bound,profile);

    \\ Checking parameters
    if(type(f) != "t_SER" && type(f) != "t_CLOSURE",error("*** Invalid type for f: ",type(f),". Has to be a modular form.***"));

    if(k_range[1] < 4,error("*** k_range[1] has to be >= 4.***"));

    if(type(k_range[2]) == "t_INT" && k_range[1] > k_range[2],error("*** Invalid range for k: ",k_range,".***"));

    \\ Initializing algorithm
    if(type(match) != "t_INT", match = if(type(f) == "t_SER", [1], 250));

    pr = if(type(f) == "t_SER", poldegree(truncate(f),'q)+1, [1]);

    profile = [match];

    phiN = eulerphi(N);
    
    k = ceil((k_range[1]-w)/phiN)*phiN + w; \\ smallest number >= k_min and cong to w mod phiN
    
    \\ Start of the main routine
    while(type(k_range[2]) != "t_INT" || k <= k_range[2],
        d = floor(k/12) + (k%12 != 2); \\ = dim M_k

        if(type(f) == "t_SER" && type(match) == "t_INT" && d + match > pr,
            print("*** Impossible to match ",match," coefficients of the series for a weight <=",k,". Increase precision.***");
            return(-1);
        );

        n = d;
        m = 1;
        while(m <= #profile,
            \\ If profile[m] is not an integer, f is a series, so pr is defined 
            bound = if(type(profile[m]) == "t_INT", d + profile[m], pr);
            kill(basis);
            basis = victor_miller_basis(k,bound,N,0); \\ The basis is not reduced
            if(m == 1,
                X = if(type(f) == "t_SER", vector(d,i,polcoeff(f,i-1,'q)), vector(d,i,f(i-1)));
                for(i=1,d-1,
                    for(j=i+1,d,
                        X[j] -= polcoeff(basis[i],j-1,'q)*X[i];
                    );
                );
            );
            while(n < bound && sum(i=1,d,X[i]*polcoeff(basis[i],n,'q)) == Mod(if(type(f) == "t_SER", polcoeff(f,n,'q), f(n)),N),
                n += 1;
            );
            \\ At this point, n <= bound.
            if(n < bound,
                break();, \\ Break the inner while
                if(m == #profile, return([k,n-d]));
            );
            m += 1;
        );
        k += phiN;
    );
    print("Impossible to match series with ",match," coefficients for a weight between ",min_k," and ",max_k," modulo ",N,".");
    return(-1);
}
addhelp(find_min_k,"find_min_k(f,w,N,{match='auto},{k_range=[4,[1]]}): Returns the minimal weight k such that the modular form f of weight w is congruent to a modular form of weight k modulo N. If the optional parameter match is set to 'auto (default), the algorithm will match all the available coefficients if f is represented by a power series and will match 200 coefficients if f is represented by a closure. Setting match to an integer forces the algorithm to stop after match coefficients match. The optional parameter k_range forces the algorithm to look for k in k_range only. The value [1] in k_range represents infinity (i.e. don't stop until a weight is found).");

find_seq(f,w,p,M,match='auto,known_seq=[]) =
{
    my(seq,k);
    if(M<=0,error("***M must be strictly positive. Recieved ",M,".***"));
    if(#known_seq == 0,
        seq = [find_min_k(f,w,p,match)];
        print("k_{",p,"^",1,"}=",seq[1]);,
        seq = known_seq;
    );
    for(n=#seq+1,M,
        k = find_min_k(f,w,p^n,match,[seq[n-1][1],[1]]);
        if(k[1] != -1,
            seq = concat(seq,[k]);
            print("k_{",p,"^",n,"}=",seq[n]);,
            print("***The sequence could not be computed up to ",p,"^",M,". ***");
            print("   but only up to ",p,"^",n-1,". ***");
            return(seq);
        );
    );
    return(seq);
}
addhelp(find_seq,"find_seq(f,w,p,M,{match='auto},{known_seq=[]}): Returns the sequence of minimal weights (see find_min_k) for the modular form f of weight w modulo p, p^2, ..., p^M. The optional parameter match has the same effect as in the find_min_k function. If known_seq is not the empty list, the algorithm will not recompute the known values.");
