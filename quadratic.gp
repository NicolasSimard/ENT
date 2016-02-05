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
