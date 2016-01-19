/* Compute the Petersson norm of delta5 using a special value of the
Symmetric square L-function attached to it. Note that delta5 is a newform. The
coefficients were computed using sage.*/

read("../computel");

w = 4; \\ Weight of delta
N = 5;

\\ initialize L-function parameters
\\ L^*(s) = A^s (gamma factors) L(s), where A = sqrt(conductor/Pi^d), where d
\\ is the number of gamma factors.
conductor = 1;              \\ exponential factor
gammaV    = [0,1,2-w];      \\ list of gamma-factors
weight    = 2*w-1;          \\ L(s)=sgn*L(weight-s)
sgn       = 1;              \\ sign in the functional equation


/* This file was produced using sage as follows:

sage: L = Newforms(Gamma0(5),4).qexp(5000).coefficients();
sage: string = "a=" + str(L) + ";"
sage: import json
sage: json.dump(string,open('coeffs.data','w'))

and then the "" were removed manually.
*/
read("coeffs.data");

print(length(a));

A = direuler(p=2,cflength(),1/(1-a[p^2]*X+a[p^2]*p^(w-1)*X^2-p^(3*w-3)*X^3));

initLdata("A[k]");        \\ L-series coefficients A(k)

print("Verifying functional equation. Error: ",errprint(checkfeq()));
print("The norm should be   = ",((4*Pi)^w/factorial(w-1)*Pi/N^2)^-1*L(w));

