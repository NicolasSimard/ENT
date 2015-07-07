/*Computation of the norm of the Delta modular form.                        */
/*The script theta_det_def.gp is too slow, so we use theorem 5.1 of Cohen's */
/*paper. Checked and it works.                                              */

k = 12

\\ Defines tau(n)
read("delta_coeff.gp");

print("Computing q-expansion...")
qExp = [];
for(i = 1, 100, qExp = concat(qExp,[tau(i)]));

f(z) = suminf(n=1, qExp[n]*exp(2*Pi*I*n*z));

\\ Define infinity and rho
oo = [1];
rho = exp(2*Pi*I/3);

V(z) = {
    local(v);
    v = [];
    for(i=0, k-2, v = concat(v,[z^i]));
    v
};

\\ We compute the integrals that go into the formula

R1 = I*intnum(t=1,[oo,2*Pi],V(I*t)*f(I*t))

R2 = intnum(t=0,1,V(t+rho)*f(t+rho))

Peter_inner = sum(n = 0, k-2,\
              (-1)^n*binomial(k-2,n)*R1[k-2-n+1]*conj(R2[n+1]))/(2*I)^(k-1)

print("computing...")
print("The answer is :",Peter_inner);
