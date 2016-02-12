isdisc(D) = return(D%4 <= 1);

control(D) = if(D%4 > 1 || D>=0,error(D," is not a negative discriminant."));

/*Every discriminant uniquely determines an order in a quadratc field. This
function returns its conductor (the index of this order in the maximal order).
*/
conductor(D) = {
    local(L,p);
    control(D);
    L=factor(D);
    p=prod(n=1,length(L~),L[,1][n]^floor(L[,2][n]/2)); \\ D=p^2m, m sq-free
    if(isdisc(D/p^2),return(p),return(p/2));
}

isprimitive(D) = return(conductor(D) == 1);

right_act(f,M) = {
    return([f[1]*M[1,1]^2+f[2]*M[1,1]*M[2,1]+f[3]*M[2,1]^2,
            2*f[1]*M[1,1]*M[1,2]+f[2]*(M[1,1]*M[2,2]+M[1,2]*M[2,1])+2*f[3]*M[2,1]*M[2,2],
            f[1]*M[1,2]^2+f[2]*M[1,2]*M[2,2]+f[3]*M[2,2]^2]
    );
}

/* Lift a matrix M in SL2(Z/NZ) to a matrix in SL2(Z).*/
liftSL2Z(M,N) = {
    local(MZ = lift(M));
    if(gcd(gcd(MZ[2,1],MZ[2,2]),N) != 1,error(M," does not belong to SL2(Z/"),N,"Z).");
    if(N<=1,return(MZ));
    local(mya,myb,myc=MZ[2,1],myd=MZ[2,2],myu=1,myv=1,myt,myg,myp);
    myg=gcd(myc,myd);
    if(myg != 1,
        myp = factor(myg)[,1];
        if(myc != 0,
            for(n=1,length(myp~),if(myc%myp[n] == 0 && myp[n] != -1,myu*=myp[n],myv*=myp[n]));
            myt=lift(chinese(Mod(1,myu),Mod(0,myv)));
            myd+=myt*N;,
            for(n=1,length(myp~),if(myd%myp[n] == 0 && myp[n] != -1,myu*=myp[n],myv*=myp[n]));
            myt=lift(chinese(Mod(1,myu),Mod(0,myv)));
            myc+=myt*N;
        );
    );
    \\ We now have a matrix mod N with coprime c and d (in Z).
    [myc,myd]=lift([myc,myd]);
    [mya,myb]=matsolvemod([myd,-myc;1,0;0,1],[0,N,N]~,[1,MZ[1,1],MZ[1,2]]~)~;
    return([mya,myb;myc,myd]);
}

 Heegner_form(f,beta,N,m1=1,m2=1) = {
    local(m,a,b,c,D,M,r,s,t,u,T);
    [a,b,c]=f;
    D=b^2-4*a*c;
    m=gcd(gcd(N,beta),(beta^2-D)/4/N);
    if(m != 1, error("Function not implemented for m>1. (D,beta,N)=(",D,",",beta,",",N,")"));
    M=[a,(b-beta)/2;(b+beta)/2,c];
    sol=matsolvemod(M,N,0,1)[2];
    if(sol[,1] != [0,0]~, [r,t]=sol[,1]~, [r,t]=sol[,2]~);
    [u,s]=matsolvemod(Mat(lift([r,-t])),N,1)~;
    T=liftSL2Z([r,s;t,u],N);
    print(T);
    return(right_act(f,T));
 }

reduced_forms(D) = {
    local(b0 = D%2, fv,a,c,zv,n);

    if(D >= 0 || D%4 > 1, return([]));

    fv = [[1,b0,(b0^2-D)/4]];
    forstep(b = b0, floor(sqrt(-D/3)), 2,
        zv=divisors((b^2-D)/4);
        n=length(zv);
        forstep(j=(n+1)\2, 2, -1,
            a=zv[j]; if (a < b, break);
            c=zv[n-j+1];
            fv = concat(fv,[[a,b,c]]);
            if(b && a != b && a != c, fv = concat(fv,[[a,-b,c]]));
        );
    );
    return(fv);
}

primitive_reduced_forms(D) = {
    local(forms, prim);
    prim = [];
    forms = reduced_forms(D);
    for(i=1, length(forms),
        if(gcd(gcd(forms[i][1],forms[i][2]),forms[i][3]) == 1,
           prim = concat(prim,[forms[i]])
        );
    );
    return(prim);
}

tau(f) = return((-f[2]+sqrt(f[2]^2-4*f[1]*f[3]))/2/f[1]);

reduced_roots(D) = {
    local(forms,taus);
    taus = [];
    forms = reduced_forms(D);
    for(i=1, length(forms),taus = concat(taus,[tau(forms[i])]));
    return(taus);
}

primitive_reduced_roots(D) = {
    local(forms,taus);
    taus = [];
    forms = primitive_reduced_forms(D);
    for(i=1, length(forms),taus = concat(taus,[tau(forms[i])]));
    return(taus);
}

class_nbr(D) = length(primitive_reduced_forms(D));

wD(D) = {
    if(D%4 > 2, error("Not a discriminant."));
    if(D == -3, return(6),
    if(D == -4, return(4),
                return(2)))
}

wQ(f) = {
    if(f[1] == f[2] && f[2] == f[3], return(6),
    if(f[1] == f[3] && f[2] == 0,    return(4),
                                     return(2)))

}

norm_class_nbr(D) = 2*class_nbr(D)/wD(D);

Hurwitz_class_nbr(D) = {
    error("Not implemented yet.")
}

genus_nbr(D) = {
    local(r,n,mu);
    if(D>=0 || D%4 > 1, return(-1));
    if(D%2 == 0, r=omega(-D)-1,r=omega(-D));
    if(D%4 == 1, return(2^(r-1)));
    n = floor(-D/4);
    if(n%4 == 3,mu=r);
    if(n%4 == 1 || n%4 == 2 || n%8 == 4,mu=r+1);
    if(n%8 == 0,mu=r+2);
    return(2^(mu-1));
}

two_torsion(D) = {
    local(forms);
    forms = primitive_reduced_forms(D);
    fv = [forms[1]];
    for(i=2, length(forms),
        if(forms[i][2] == 0 || forms[i][1] == forms[i][3] || abs(forms[i][2]) == forms[i][1],
        fv = concat(fv,[forms[i]]);
        );
    );
    return(fv);
}

\\ Returns the Chowla-selberg period of discriminant D, as defined in 1-2-3 of
\\ modular forms by Zagier. Tested for D = -4.
CSperiod(D) = {
    if(D%4 > 2, error("Not a discriminant."));
    return(prod(j=1,abs(D)-1,gamma(j/abs(D))^kronecker(D,j))^(wD(D)/4/class_nbr(D))/sqrt(2*Pi*abs(D)));
}
