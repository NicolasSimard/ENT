/*Just a small script containing the function that I use very often.
*/

/* To define the Eisenstein series, we follow Zagier's convention:

G_k(s) = 1/2*sum_{(m,n)\neq(0,0)}(mz+n)^-k

E_k(s)=G_k(s)/zeta(k)

GG_k(s) = (k-1)!/(2*Pi*I)^kG_k(s) (all non-constant Fourrier coeff are int)
*/
GG(k,s) = -bernfrac(k)/2/k+suminf(n=1,sigma(n,k-1)*exp(2*Pi*I*n*s));

GG_qexp(k,prec) = -bernfrac(k)/2/k+sum(n=1,prec,sigma(n,k-1)*q^n)+O(q^(prec+1));

E(k,s) = -2*k/bernfrac(k)*GG(k,s);

E_qexp(k,prec) = -2*k/bernfrac(k)*GG_qexp(k,prec);

/* Note: recall that the coefficients of delta are given by tau(n), where
? tau(n) = tau(n) = (5*sigma(n,3)+7*sigma(n,5))*n/12-35*sum(k=1,n-1,(6*k-4*(n-k))*sigma(k,3)*sigma(n-k,5));

However, it is much faster to do

? D = delta_qexp(1000);
? L = vector(1000,n,polcoeff(truncate(D),n,q));

than to do

? L = vector(1000,n,tau(n));
*/
delta_qexp(prec) = (E_qexp(4,prec)^3-E_qexp(6,prec)^2)/1728;

/* Takes 355 ms to compute 1000 coefficients.*/
j_qexp(prec) = E_qexp(4,prec+2)^3/delta_qexp(prec+2);

qexp_coeff(f) = {
    local(v,p);
    v = valuation(f,q);
    p = truncate(q^(1-v)*f);
    return(vector(poldegree(p),n,polcoeff(p,n,q)));
    return(vector)
}
