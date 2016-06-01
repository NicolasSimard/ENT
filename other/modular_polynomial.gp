\r ../Modform.gp

j_gamma(a,b,d,verbose=0) = {
    my(j,p,pr=default(seriesprecision));
    if(verbose,print("(",a,",",b,",",d"):"));
    j = ellj('q);
    p = lift(Mod(sum(n=-1,pr-2,polcoeff(j,n,'q)*'zd^((b*n)%d)*'qd^(a*n)),polcyclo(d,'zd)));
    return(p + O('qd^(a*pr+1)));
}

inner_pol(a,d,verbose=0) = {
    my(p,p_norm);
    p = prod(b=0,d-1,'X-j_gamma(a,b,d,verbose));
    p = lift(Mod(p,polcyclo(d,'zd)));
    p_norm=substpol(p,'qd^d,'q); \\ Replaces qd^d by q in p.
    return(p_norm);
}

/* Returns the modular polynomial psi_m(X,Y).Note that the precision \p and \ps
have to be high enough since we are rounding and using jpol.*/
psi_m(m,verbose=0) = {
    my(p=1,deg);
    fordiv(m,d,p*=inner_pol(m/d,d,verbose));
    deg=poldegree(truncate(p),'X);
    return(Pol(vector(deg+1,n,jpol(polcoeff(p,deg+1-n,'X)+O('q))),'Y));
}

/*Returns the modular polynomial psi_m(X,X). Note that it has degree
sigma^+(m)=sum_{d|m} max(d,m/d). */
psi_mX(m,verbose=0) = subst(psi_m(m,verbose),'Y,'X);

H_d(d) =
{
    if(issquare(d) || issquare(d/3), error("Not implemented yet!"));
    \\print("Computing fd.");
    my(Hd = class_nbr(d), mffd = fd(d), p);
    \\print("Computing Borcherds product.");
    p = 'q^-Hd*prod(n=1,Hd,(1-'q^n+O('q^(Hd+1)))^polcoeff(mffd,n^2,'q));
    \\print("Computing jpol.");
    jpol(p);
}
