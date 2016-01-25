/*Computing the Petersson norm using the special values of the symmetric
L-function of Delta.

Tested up to the number of decimals in Cohen's paper.
*/

read("../computel");           \\ read the ComputeL package

w = 12; \\ Weight of delta

\\ initialize L-function parameters
\\ L^*(s) = A^s (gamma factors) L(s), where A = sqrt(conductor/Pi^d), where d
\\ is the number of gamma factors.
conductor = 1;              \\ exponential factor
gammaV    = [0,1,2-w];      \\ list of gamma-factors
weight    = 2*w-1;          \\ L(s)=sgn*L(weight-s)
sgn       = 1;              \\ sign in the functional equation

read("delta_coeff.gp");     \\define tau(k)
A(k) = sumdiv(k, d, (-1)^bigomega(d)*d^(w-1)*tau(k/d)^2);

initLdata("A(k)");        \\ L-series coefficients A(k)

print("Verifying functional equation. Error: ",errprint(checkfeq()));
print("The Petersson norm is: ",2/Pi*factorial(w-1)/(4*Pi)^w*L(w));
