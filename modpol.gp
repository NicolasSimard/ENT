read("modforms.gp");

j_gamma(a,b,d,prec) = {
    local(j);
    j = j_qexp(prec);
    return(sum(n=-1,prec,polcoeff(j,n,q)*exp(2*Pi*I*b*n/d)*qd^(a*n))+O(qd^(a*prec+1)));
}

inner_pol(a,d,prec) = {
    local(p,p_norm);
    p = round(real(prod(b=0,d-1,X-j_gamma(a,b,d,prec))));
    lowest =valuation(p,qd); \\ Is a multiple of d
    highest=a*prec;
print("lowest:",lowest,"    highest:",highest);
    p_norm = sum(n=lowest/d,floor(highest/d),polcoeff(p,n*d,qd)*q^n)+O(q^(highest+1));
\\print(p);
\\print(p_norm);
    return(p_norm);
}

psi_Xj(m) = {
    local(p=1);
    fordiv(m,d,p*=inner_pol(m/d,d,100));
    return(p);
}

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
