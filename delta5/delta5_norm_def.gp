/***Computing the Petersson norm of delta5 directly from the definition.  ***/
/***This code gives the right answer!                                     ***/

delta5(z) = (eta(z,1)*eta(5*z,1))^4; \\The 1 indicates that we want the "true" eta (with q^(1/24))

oo = [1];                            \\We define infinity

\\ We split the integral over a fundamental domain for Gamma_0(5) into
\\ integrals over a fundamental domain for SL_2(Z) by using coset reps and
\\ then summing.

k = 4
p = 5

sub(i) = intnum(x = -1/2, 1/2,intnum(y = (1-x^2)^(1/2),[oo,4*Pi],norm(delta5((x+y*I)/(i*(x+y*I)+1))*(i*(x+I*y)+1)^-k)*y^(k-2)));

sub0 = intnum(x = -1/2, 1/2,intnum(y = (1-x^2)^(1/2),[oo,4*Pi],norm(delta5(-1/(x+y*I))*(x+I*y)^-k)*y^(k-2)));

print("Computing the norm using the definition.");
N = (sub0 + sum(i = 0, p-1, sub(i)))/(p+1);
print("The norm is: ", N)