/* Compute the Petersson norm of delta5 using a special value of the
Symmetric square L-function attached to it. Note that delta5 is a newform. The
coefficients were computed using sage.*/

read("../computel");

w = 4; \\ Weight of delta
N = 5;

\\ initialize L-function parameters
\\ L^*(s) = A^s (gamma factors) L(s), where A = sqrt(conductor/Pi^d), where d
\\ is the number of gamma factors.
conductor = N;              \\ exponential factor
gammaV    = [0,1,2-w];      \\ list of gamma-factors
weight    = 2*w-1;          \\ L(s)=sgn*L(weight-s)

/* This file was produced using sage as follows:
sage
: L = Newforms(Gamma0(5),4).qexp(5000).coefficients();
sage: string = "a=" + str(L) + ";"
sage: import json
sage: json.dump(string,open('coeffs.data','w'))

and then the "" were removed manually.
*/
read("partial_coeff.data");

print(cflength());

A = direuler(p=2,cflength(),1/(1-a[p^2]*X+a[p^2]*p^(w-1)*X^2-p^(3*w-3)*X^3));

initLdata("A[k]");        \\ L-series coefficients A(k)


sgn       = X;              \\ sign in the functional equation
sgneq     = Vec(checkfeq(1.1));       \\ checkfeq() returns c1+X*c2, should be 0
sgn       = -sgneq[2]/sgneq[1];    \\ hence solve for the sign

print("Verifying functional equation. Error: ",errprint(checkfeq()));

index = N;
forprime(p=2,N,if(N%p == 0,index = index*(1+p^-1)));

norme = ((4*Pi)^w/factorial(w-1)*Pi/N^2/eulerphi(N))^-1/index*L(w);
print("The norm should be   = ",norme);

