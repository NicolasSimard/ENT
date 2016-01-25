read("../computel");

w = 12; \\ Weight of delta
N = 1;

\\ initialize L-function parameters
\\ L^*(s) = A^s (gamma factors) L(s), where A = sqrt(conductor/Pi^d), where d
\\ is the number of gamma factors.
conductor = 1;              \\ exponential factor
gammaV    = [0,1,2-w];      \\ list of gamma-factors
weight    = 2*w-1;          \\ L(s)=sgn*L(weight-s)
sgn       = 1;              \\ sign in the functional equation

read("delta_coeff.gp");

A = direuler(p=2,cflength(),1/(1-tau(p^2)*X+tau(p^2)*p^(w-1)*X^2-p^(3*w-3)*X^3));

initLdata("A[k]");        \\ L-series coefficients A(k)

print("Verifying functional equation. Error: ",errprint(checkfeq()));
print("The norm should be   = ",((4*Pi)^w/factorial(w-1)*Pi/2/N^2)^-1*L(w));

