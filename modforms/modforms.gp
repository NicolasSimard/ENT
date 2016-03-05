/*Just a small script to compute a few q-expansions.*/

/* To define the Eisenstein series, we follow Zagier's convention:

G_k(s) = 1/2*sum_{(m,n)\neq(0,0)}(mz+n)^-k

E_k(s)=G_k(s)/zeta(k) (leading coefficient is 1)

GG_k(s) = (k-1)!/(2*Pi*I)^kG_k(s) (all non-constant Fourrier coeff are int)
*/
GG(k:small,s) = -bernfrac(k)/2/k+suminf(n=1,sigma(n,k-1)*exp(2*Pi*I*n*s));

GG_qexp(k:small,pr:small) = -bernfrac(k)/2/k+Ser(concat([0],vector(pr-1,n,sigma(n,k-1))),'q);

E(k:small,s) = -2*k/bernfrac(k)*GG(k,s);

E_qexp(k:small,pr:small) = -2*k/bernfrac(k)*GG_qexp(k,pr);

/* Auxilary function to compute the q-expansion of Delta.*/
eta3_qexp(pr:small) = 
{
    my(d);
    Ser(vector(pr,n,if(issquare(1+8*(n-1),&d),(-1)^((-1+d)/2)*d,0)),'q);
}

/* This function is much faster than this one
theta_qexp2(pr) = 1+sum(n=1,floor(sqrt(pr)),2*q^(n^2)) + O(q^pr); 
and uses much less memory!
*/
theta_qexp(pr:small) = Ser(vector(pr,n,2*issquare(n-1)),'q)-1;

/* Note: recall that the coefficients of delta are given by tau(n), where
? tau(n) = tau(n) = (5*sigma(n,3)+7*sigma(n,5))*n/12-35*sum(k=1,n-1,(6*k-4*(n-k))*sigma(k,3)*sigma(n-k,5));

However, it is much faster to do

? D = delta_qexp(1000);
? L = vector(1000,n,polcoeff(truncate(D),n,q));

than to do

? L = vector(1000,n,tau(n));
*/
delta_qexp(pr:small) = 'q*sqr(sqr(sqr(eta3_qexp(pr-1))));

/* Takes 213 ms to compute 1000 coefficients.*/
j_qexp(pr:small) = E_qexp(4,pr+2)^3/delta_qexp(pr+2);
