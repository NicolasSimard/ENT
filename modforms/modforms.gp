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

/* Auxilary function to compute the q-expansion of Delta.*/
theta1(prec) = sum(n=0,floor((-1+sqrt(1+8*prec))/2),(-1)^n*(2*n+1)*q^(n*(n+1)/2))+O(q^(prec+1));

/* Note: recall that the coefficients of delta are given by tau(n), where
? tau(n) = tau(n) = (5*sigma(n,3)+7*sigma(n,5))*n/12-35*sum(k=1,n-1,(6*k-4*(n-k))*sigma(k,3)*sigma(n-k,5));

However, it is much faster to do

? D = delta_qexp(1000);
? L = vector(1000,n,polcoeff(truncate(D),n,q));

than to do

? L = vector(1000,n,tau(n));
*/
delta_qexp(prec) = q*theta1(prec-1)^8;

/* The traditionnal definition of Delta. Much slower.*/
delta_qexp2(prec) = (E_qexp(4,prec)^3-E_qexp(6,prec)^2)/1728;

/* Takes 355 ms to compute 1000 coefficients.*/
j_qexp(prec) = E_qexp(4,prec+2)^3/delta_qexp(prec+2);

qexp_coeff(f) = {
    local(v,p);
    v = valuation(f,q);
    p = truncate(q^(1-v)*f);
    return(vector(poldegree(p),n,polcoeff(p,n,q)));
}

victor_miller_basis(k,prec=10,N=0,notred=0) = {
    local(n,e,ls,A,E6,E6_squared,D,Eprod,Dprod);
    if(k%2 == 1,error("The weight ",k," must be even."));

    e=k%12 + 12*(k%12 == 2);
    n=(k-e)/12;

    E6 = E_qexp(6,prec);
    A = if(e==0,1,E_qexp(e,prec)); \\ Slight difference with e=6
    if(polcoeff(A,0,q) == -1, A=-A);
    D = delta_qexp(prec);

    if(N>0,E6 = Mod(E6,N); A = Mod(A,N); D = Mod(D,N));

    E6_squared = E6^2;
    Eprod = E6_squared;
    Dprod = D;

    ls = vector(n+1,i,A);
    for(i=1,n,
        ls[n-i+1] *= Eprod;
        ls[i+1] *= Dprod;
        Eprod *= E6_squared;
        Dprod *= D;
    );
    \\ At this point, the basis is upper triangular
    if(notred,return(ls));
    for(i=1,n,
        for(j=1,i,
            ls[j] = ls[j] - polcoeff(ls[j],i,q)*ls[i+1];
        );
    );
    return(ls);
}

