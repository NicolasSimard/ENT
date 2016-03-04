\r modforms.gp

victor_miller_basis(k,pr=10,N=0,reduce=1) = 
{
    my(n,e,ls,A,E6,E6_squared,D,Eprod,Dprod);
    if(k%2 == 1,error("The weight ",k," must be even."));

    e=k%12 + 12*(k%12 == 2);
    n=(k-e)/12;

    if(pr <= n,
        ls = vector(n+1);
        ls[1] = 1+ O(q^pr);
        for(i=2,pr,ls[i] = q^i + O(q^pr));
        for(i=pr+1,n+1,ls[i] = O(q^pr));
        return(ls);
    );

    E6 = E_qexp(6,pr);
    A = if(e==0,1,E_qexp(e,pr)); \\ Slight difference with e=6
    if(polcoeff(A,0,q) == -1, A=-A);
    D = delta_qexp(pr);

    if(N>0,E6 = E6*Mod(1,N); A = A*Mod(1,N); D = D*Mod(1,N));
    \\if(N>0,E6 = Mod(E6,N); A = Mod(A,N); D = Mod(D,N));

    E6_squared = sqr(E6) + O(q^pr);
    Eprod = E6_squared;
    Dprod = D;

    ls = vector(n+1,i,A);
    for(i=1,n,
        ls[n-i+1] = ls[n-i+1]*Eprod + O(q^pr);
        ls[i+1] = ls[i+1]*Dprod + O(q^pr);
        Eprod = Eprod*E6_squared + O(q^pr);
        Dprod = Dprod*D + O(q^pr);
    );
    \\At this point the basis is upper triangular
    if(reduce == 0,return(ls));
    for(i=1,n,
        for(j=1,i,
            ls[j] = ls[j] - polcoeff(ls[j],i,q)*ls[i+1];
        );
    );
    return(ls);
}

