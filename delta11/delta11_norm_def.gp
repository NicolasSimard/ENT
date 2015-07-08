/* Compute the norm of delta11 using the definition of the Petersson inner  */
/* product. \p 25 gives 6 decimals.                                         */

\\ Define the parameters
p = 11;
k = 2;

\\ Define the actual function
f(z) = (eta(z,1)*eta(11*z,1))^2; \\ f === delta11

oo = [1]

sub(i) = intnum(x = -1/2, 1/2,intnum(y = (1-x^2)^(1/2),[oo,4*Pi],norm(f((x+y*I)/(i*(x+y*I)+1))*(i*(x+I*y)+1)^-k)*y^(k-2)));

sub0   = intnum(x = -1/2, 1/2,intnum(y = (1-x^2)^(1/2),[oo,4*Pi],norm(f(-1/(x+y*I))*(x+I*y)^-k)*y^(k-2)));

print("Computing the norm using the definition.");
N = (sub0 + sum(i = 0, p-1, sub(i)))/(p+1);
print("The norm is: ", N)