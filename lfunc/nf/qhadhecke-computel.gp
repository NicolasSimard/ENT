/* Define the Hecke L-function attached to Hecke character of an imaginary
quadratic field. The Hecke character must have conductor 1 for the moment.

Tested with the values in Watkins's paper on computations with Hecke chars.
*/

read("../computel");
read("Qhc.gp");

ell=50;
K = bnfinit('x^2+23);
T = [4*ell,0];
comp = [1];

qhcdata = qhcinit(K);

\\ Data of the L-function
conductor   = abs(K.disc);
gammaV      = [0,1];
weight      = T[1]+T[2]+1;
sgn         = 'X; \\ Sign will be computed automatically
Lpoles      = [];
Lresidues   = [];

\\Computing the coefficients of the Dirichlet series
print("Computing coefficients");
idealsofnorm = ideallist(K,cflength());
a = vector(cflength(),k,sum(i=1,#idealsofnorm[k],qhceval(qhcdata,[comp,T],idealsofnorm[k][i])));

\\ Initialise the L-function
print("Initializing");
initLdata("a[k]",,"conj(a[k])");

\\ Compute the sign
print("computing sign");
sgneq = Vec(checkfeq());       \\ checkfeq() returns c1+X*c2, should be 0
sgn   = -sgneq[2]/sgneq[1];    \\ hence solve for the sign