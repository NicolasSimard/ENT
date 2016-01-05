read("computel");           \\ read the ComputeL package
read("delta_coeff.gp");     \\define tau(k)

coefgrow(n) = 2*n^(11/2);

\\ initialize L-function parameters
conductor = 1;              \\ exponential factor
gammaV    = [0,1];          \\ list of gamma-factors
weight    = 12;             \\ L(s)=sgn*L(weight-s)
sgn       = 1;              \\ sign in the functional equation

initLdata("tau(k)");        \\ L-series coefficients a(k)

print("Verifying functional equation. Error: ",errprint(checkfeq()))