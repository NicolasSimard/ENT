/* Compute the norm of delta11 using the definition of the Petersson inner  */
/* product. \p 25 gives 6 decimals.                                         */

\\ Define the parameters
p = 11;
k = 2;

\\ Define the actual function
f(z) = (eta(z,1)*eta(11*z,1))^2; \\ f === delta11

\\ Deffine the coefficients of delta11. Obtained from sage via
\\ CuspForms(11,2).q_expansion_basis(100)[0].coefficients(), since
\\ the space S_2(\Gamma_0(11)) has dimension 1.
read("delta11_coeffs.gp"); \\ This defines the array coeff

print(length(coeff));

ff(z) = suminf(n=1, coeff[n]*exp(2*Pi*I*n*z))

oo = [1]

sub(i) = intnum(x = -1/2, 1/2,intnum(y = (1-x^2)^(1/2),[oo,4*Pi],\
          norm(ff((x+y*I)/(i*(x+y*I)+1))*(i*(x+I*y)+1)^-k)*y^(k-2)));

sub0   = intnum(x = -1/2, 1/2,intnum(y = (1-x^2)^(1/2),[oo,4*Pi],\
          norm(ff(-1/(x+y*I))*(x+I*y)^-k)*y^(k-2)));

print("Computing the norm using the definition.");
N = (sub0 + sum(i = 0, p-1, sub(i)))/(p+1);
print("The norm is: ", N)