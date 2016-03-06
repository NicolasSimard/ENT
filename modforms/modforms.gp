/* This script defines a few modular forms. These modular forms are represented
as functions sending integers n to a_n(f), the nth Fourier  coefficient of f.

Some modular forms can also be represented as Fourrier series directly. This is
often faster than computing every Fourrier coefficient one after the other. For
example, delta_qexp(1000) is much faster than Ser(vector(1000,n,delta(n))).*/

/* To define the Eisenstein series, we follow Zagier's convention:

G_k(s) = 1/2*sum_{(m,n)\neq(0,0)}(mz+n)^-k

E_k(s)=G_k(s)/zeta(k) (leading coefficient is 1)

GG_k(s) = (k-1)!/(2*Pi*I)^kG_k(s) (all non-constant Fourrier coeff are int)
*/

\\GG(k:small,s) = -bernfrac(k)/2/k+suminf(n=1,sigma(n,k-1)*exp(2*Pi*I*n*s));

GG(k:small) = 
{
    n -> if(n == 0, -bernfrac(k)/2/k, sigma(n,k-1));
}
addhelp(GG,"GG(k): Returns Eisenstein series GG of weight k. Modular forms are\nrepresented as closures, sending n to the nth Fourrier coefficient.");

GG_qexp(k:small,pr:small) = 
{
    -bernfrac(k)/2/k+Ser(concat([0],vector(pr-1,n,sigma(n,k-1))),'q);
}
addhelp(GG_qexp,"GG_qexp(k,pr): Returns the q-expansion of the Eisenstein series\nGG of weight k.");

\\E(k:small,s) = -2*k/bernfrac(k)*GG(k,s);

E(k:small) =
{
    n -> if(n == 0, 1, -2*k/bernfrac(k)*sigma(n,k-1));
}
addhelp(E,"E(k): Returns Eisenstein series E of weight k. Modular forms are\nrepresented as closures, sending n to their nth Fourrier coefficient.");

E_qexp(k:small,pr:small) =
{
    -2*k/bernfrac(k)*GG_qexp(k,pr);
}
addhelp(E_qexp,"E_qexp(k,pr): Returns the q-expansion of the Eisenstein series\nE of weight k.");

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
theta_qexp(pr:small) =
{
    Ser(vector(pr,n,2*issquare(n-1)),'q)-1;
}
addhelp(theta_qexp,"theta_qexp(k): Returns the q-expansion of the theta series of weight 1/2\n up to precision pr.");

delta(n) = polcoeff(delta_qexp(n+1),n,'q);
addhelp(delta,"delta: Returns the modular form delta. Modular forms are represented\nas functions sending n to their nth Fourrier coefficient.");


delta_qexp(pr:small) =
{
    'q*sqr(sqr(sqr(eta3_qexp(pr-1))));
}
addhelp(delta_qexp,"delta_qexp(k): Returns the q-expansion of the delta series up to precision pr.");

/* Takes 213 ms to compute 1000 coefficients.*/
j_qexp(pr:small) =
{
    E_qexp(4,pr+2)^3/delta_qexp(pr+2);
}
addhelp(j_qexp,"j_qexp(pr): Returns the q-expansion of the j-invariant up to precision pr.");

clos2qexp(f,pr:small) = 
{
    return(Ser(vector(pr,n,f(n-1)),'q));
}
addhelp(clos2qexp,"cols2qexp(f,pr): Return the q-expansion attached to the closure representation of f up to precision pr.");

qexp2clos(f) =
{
    n -> polcoeff(f,n,'q);
}
addhelp(qexp2clos,"qexp2clos(f): Returns the closure attached to the q-expansion of f.");
