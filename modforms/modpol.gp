read("modforms.gp");

j_gamma(a,b,d,prec) = {
    local(j,p);
    j = j_qexp(prec);
    p = lift(Mod(sum(n=-1,prec,polcoeff(j,n,q)*zd^((b*n)%d)*qd^(a*n)),polcyclo(d,zd)));
    return(p + O(qd^(a*prec+1)));
}

inner_pol(a,d,prec) = {
    local(p,p_norm);
    p = prod(b=0,d-1,X-j_gamma(a,b,d,prec));
    p = lift(Mod(p,polcyclo(d,zd)));
    p_norm=substpol(p,qd^d,q); \\ Replaces qd^d by q in p.
    return(p_norm);
}

/* Returns the modular polynomial psi_m(X,Y). The parameter prec determines
the number of coefficients of the q-expansion of j that will be computed.
If too low, it becomes impossible to transform the coefficients of psi_m(X,j)
into polynomials in j. Note also that the precision has the be high enough
since we are rounding.*/
psi_m(m,prec=200) = {
    local(p=1,deg);
    fordiv(m,d,p*=inner_pol(m/d,d,prec));
    deg=poldegree(p,X);
    return(Pol(vector(deg+1,n,j_pol(polcoeff(p,deg+1-n,X))),X));
}

/*Returns the modular polynomial psi_m(X,X). Note that it has degree
sigma^+(m)=sum_{d|m} max(d,m/d). */
psi_mX(m,prec=200) = subst(psi_m(m,prec),Y,X);

/* Given a q-expansion f(q) with enough coefficients, returns a polynomial
g(Y) such that g(j)=f, where j is the j-function.
*/
j_pol(f) = {
    local(k=-valuation(f,q));
    j=j_qexp(k-1);
    js=vector(k+1,n,j^(k+1-n));
    M=matrix(k+1,k+1,n,m,polcoeff(js[m],n-k-1,q));
    B=vector(k+1,n,polcoeff(f,n-k-1,q));
    return(Pol(matsolve(M,B~),Y));
}
