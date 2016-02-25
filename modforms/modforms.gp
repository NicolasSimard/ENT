/*Just a small script to compute a few q-expansions.*/

/* To define the Eisenstein series, we follow Zagier's convention:

G_k(s) = 1/2*sum_{(m,n)\neq(0,0)}(mz+n)^-k

E_k(s)=G_k(s)/zeta(k) (leading coefficient is 1)

GG_k(s) = (k-1)!/(2*Pi*I)^kG_k(s) (all non-constant Fourrier coeff are int)
*/
GG(k,s) = -bernfrac(k)/2/k+suminf(n=1,sigma(n,k-1)*exp(2*Pi*I*n*s));

GG_qexp(k,prec) = -bernfrac(k)/2/k+Ser(concat([0],vector(prec-1,n,sigma(n,k-1))),q);

E(k,s) = -2*k/bernfrac(k)*GG(k,s);

E_qexp(k,prec) = -2*k/bernfrac(k)*GG_qexp(k,prec);

/* Auxilary function to compute the q-expansion of Delta.*/
eta3_qexp(prec) = {
    local(k,d);
    Ser(vector(prec,n,if(issquare(1+8*(n-1),&d),(-1)^((-1+d)/2)*d,0)),q);
}

/* This function is much faster than this one
theta_qexp2(prec) = 1+sum(n=1,floor(sqrt(prec)),2*q^(n^2)) + O(q^prec); 
and uses much less memory! Avoid it!
*/
theta_qexp(prec) = Ser(vector(prec,n,2*issquare(n-1)),q)-1;

/* Note: recall that the coefficients of delta are given by tau(n), where
? tau(n) = tau(n) = (5*sigma(n,3)+7*sigma(n,5))*n/12-35*sum(k=1,n-1,(6*k-4*(n-k))*sigma(k,3)*sigma(n-k,5));

However, it is much faster to do

? D = delta_qexp(1000);
? L = vector(1000,n,polcoeff(truncate(D),n,q));

than to do

? L = vector(1000,n,tau(n));
*/
delta_qexp(prec) = q*sqr(sqr(sqr(eta3_qexp(prec-1))));

/* Takes 213 ms to compute 1000 coefficients.*/
j_qexp(prec) = E_qexp(4,prec+2)^3/delta_qexp(prec+2);

victor_miller_basis(k,prec=10,N=0,reduce=1) = {
    local(n,e,ls,A,E6,E6_squared,D,Eprod,Dprod);
    if(k%2 == 1,error("The weight ",k," must be even."));

    e=k%12 + 12*(k%12 == 2);
    n=(k-e)/12;

    if(prec <= n,
        ls = vector(n+1);
        ls[1] = 1+ O(q^prec);
        for(i=2,prec,ls[i] = q^i + O(q^prec));
        for(i=prec+1,n+1,ls[i] = O(q^prec));
        return(ls);
    );

    E6 = E_qexp(6,prec);
    A = if(e==0,1,E_qexp(e,prec)); \\ Slight difference with e=6
    if(polcoeff(A,0,q) == -1, A=-A);
    D = delta_qexp(prec);

    if(N>0,E6 = E6*Mod(1,N); A = A*Mod(1,N); D = D*Mod(1,N));
    \\if(N>0,E6 = Mod(E6,N); A = Mod(A,N); D = Mod(D,N));

    E6_squared = sqr(E6) + O(q^prec);
    Eprod = E6_squared;
    Dprod = D;

    ls = vector(n+1,i,A);
    for(i=1,n,
        ls[n-i+1] = ls[n-i+1]*Eprod + O(q^prec);
        ls[i+1] = ls[i+1]*Dprod + O(q^prec);
        Eprod = Eprod*E6_squared + O(q^prec);
        Dprod = Dprod*D + O(q^prec);
    );
    \\ At this point, the basis is upper triangular
    if(reduce == 0,return(ls));
    for(i=1,n,
        for(j=1,i,
            ls[j] = ls[j] - polcoeff(ls[j],i,q)*ls[i+1];
        );
    );
    return(ls);
}

