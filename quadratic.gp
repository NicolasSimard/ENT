/*1) reduced_forms(D): return a list of triples representing the reduced forms
of discriminant D.

2) theta_qexps(D, k, N): return a list of the coefficients of the theta series
attached to the quadratic field of discriminant D and with parameter k. Note
that the weight of these theta series is not k.
*/

reduced_forms(D) = {
    local(b0 = D%2, fv,a,c,zv);

    if(D >= 0 || D%4 > 1, return([]));

    fv = [[1,b0,(b0^2-D)/4]];
    forstep(b = b0, floor(sqrt(-D/3)), 2,
        zv=divisors((b^2-D)/4);
        n=length(zv);
        forstep(j=(n+1)\2, 2, -1,
            a=zv[j]; if (a < b, break);
            c=zv[n-j+1];
            if(gcd(gcd(a,b),c) != 1, next);
            fv = concat(fv,[[a,b,c]]);
            if(b && a != b && a != c, fv = concat(fv,[[a,-b,c]]));
        );
    );
    fv
}

reduced_roots(D) = {
    local(forms,taus);
    taus = [];
    forms = reduced_forms(D);
    for(i=1, length(forms),
        taus = concat(taus,[(-forms[i][2]+sqrt(D))/2/forms[i][1]]);
    );
    taus
}

class_nbr(D) = {
    length(reduced_forms(D))
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
    2^(mu-1)
}

two_torsion(D) = {
    local(forms);
    forms = reduced_forms(D);
    fv = [forms[1]];
    for(i=2, length(forms),
        if(forms[i][2] == 0 || forms[i][1] == forms[i][3] || abs(forms[i][2]) == forms[i][1],
        fv = concat(fv,[forms[i]]);
        );
    );
    fv
}

