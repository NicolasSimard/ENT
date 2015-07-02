/***Computing the Petersson norm of delta directly from the definition.   ***/
/***This code gives the right answer!                                     ***/

delta(z) = eta(z,1)^24; \\The 1 indicates that we want the "true" eta (with q^(1/24))

oo = [1];               \\We define infinity

print("Computing the norm using the definition.");
N = intnum(x = -1/2, 1/2, intnum(y = (1-x^2)^(1/2),[oo,4*Pi], norm(delta(x+y*I))*y^10));
print("The norm is: ", N);

\\ We now compute delta using its q-expansion, for fun!
read("delta_coeff.gp");

my_delta(z) = suminf(n=1, tau(n)*exp(2*Pi*I*n*z));

print("Computing the norm using the definition.");
N2 = intnum(x = -1/2, 1/2, intnum(y = (1-x^2)^(1/2),[oo,4*Pi], norm(my_delta(x+y*I))*y^10));
print("The norm is: ", N2);