\r modforms.gp

j_gamma(a,b,d,pr,verbose=0) = {
    my(j,p);
    if(verbose,print("(",a,",",b,",",d"):",pr));
    j = j_qexp(pr);
    p = lift(Mod(sum(n=-1,pr-1,polcoeff(j,n,'q)*'zd^((b*n)%d)*'qd^(a*n)),polcyclo(d,'zd)));
    return(p + O('qd^(a*pr+1)));
}

inner_pol(a,d,pr,verbose=0) = {
    my(p,p_norm);
    p = prod(b=0,d-1,'X-j_gamma(a,b,d,pr,verbose));
    p = lift(Mod(p,polcyclo(d,'zd)));
    p_norm=substpol(p,'qd^d,'q); \\ Replaces qd^d by q in p.
    return(p_norm);
}

/* Returns the modular polynomial psi_m(X,Y). The parameter pr determines
the number of coefficients of the q-expansion of j that will be computed.
If too low, it becomes impossible to transform the coefficients of psi_m(X,j)
into polynomials in j. Note also that the precision has the be high enough
since we are rounding.*/
psi_m(m,pr=200,verbose=0) = {
    my(p=1,deg);
    fordiv(m,d,p*=inner_pol(m/d,d,pr,verbose));
    deg=poldegree(truncate(p),'X);
    return(Pol(vector(deg+1,n,j_pol(polcoeff(p,deg+1-n,'X)+O('q))),'X));
}

/*Returns the modular polynomial psi_m(X,X). Note that it has degree
sigma^+(m)=sum_{d|m} max(d,m/d). */
psi_mX(m,pr=200,verbose=0) = subst(psi_m(m,pr,verbose),'Y,'X);

/* Given a q-expansion f(q) with enough coefficients, returns a polynomial
g(Y) such that g(j)=f, where j is the j-function.
*/
j_pol(f) = {
    my(k=-valuation(f,'q));
    if(k <=  0, return(f));
    j=j_qexp(k);
    js=vector(k+1,n,j^(k+1-n));
    M=matrix(k+1,k+1,n,m,polcoeff(js[m],n-k-1,'q));
    B=vector(k+1,n,polcoeff(f,n-k-1,'q));
    return(Pol(matsolve(M,B~),'Y));
}
