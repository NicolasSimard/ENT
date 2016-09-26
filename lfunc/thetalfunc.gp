/*Define the L-function of the theta series attached to the quadratic fields.

Tested for a few quadratic fields and it works.
*/

read("../computel");                 \\ read the ComputeL package
read("theta_function.gp");

\\ Define the parameters. Can be owerwritten.
p = 47;
Darmonk = 3;
f = reduced_forms(-p)[3]
verbose = 1;
                            \\ initialize L-function parameters
conductor = p;              \\ exponential factor
gammaV    = [0,1];          \\ list of gamma-factors
weight    = 2*Darmonk + 1;  \\ L(s)=sgn*L(weight-s)
sgn       = 1;              \\ sign in the functional equation

if(verbose, print("Computing q-expansion."));
coeff = theta_qexp(f, Darmonk, 2*cflength());

initLdata("coeff[k]");        \\ L-series coefficients a(k)

print("Verifying functional equation for p=",p,", k=",Darmonk,"error: ",errprint(checkfeq()));
print("L(1)       = ",lval = L(1));
print(" (check)   = ",lval2 = L(1,1.1),"  (err=",errprint(lval-lval2),")");
