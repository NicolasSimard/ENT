/* Compute the Petersson norm of delta5 using a special value of the
Symmetric square L-function attached to it. Note that delta5 is a newform. The
coefficients were computed using sage.*/

read("../computel");

k = 2; \\ Weight of delta
N = 11;

\\ initialize L-function parameters
\\ L^*(s) = A^s (gamma factors) L(s), where A = sqrt(conductor/Pi^d), where d
\\ is the number of gamma factors.
conductor = N^2;              \\ exponential factor
gammaV    = [0,1,k-1,k];      \\ list of gamma-factors
weight    = 1;              \\ L(s)=sgn*L(weight-s)
Lpoles    = [1];
sgn       = 1;

read("modified-Rankin-coeffs.data");

print(cflength());

initLdata("a[k]");        \\ L-series coefficients A(k)

/*sgn       = X;              \\ sign in the functional equation
sgneq     = Vec(checkfeq(1.1));       \\ checkfeq() returns c1+X*c2, should be 0
sgn       = -sgneq[2]/sgneq[1];    \\ hence solve for the sign
*/

print("Verifying functional equation. Error: ",errprint(checkfeq()));

/*index = N;
forprime(p=2,N,if(N%p == 0,index = index*(1+p^-1)));

norme = ((4*Pi)^w/factorial(w-1)*Pi/N^2/eulerphi(N))^-1/index*L(w);
other_norm = 2/Pi*factorial(w-1)*N/(4*Pi)^w*index*L(1);
print("The norm should be   = ",norme);
print("or                   = ",other_norm);*/
