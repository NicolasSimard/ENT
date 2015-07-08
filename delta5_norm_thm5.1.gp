/*Computation of the norm of delta5.                                        */
/*The script theta_det_def.gp is too slow, so we use theorem 5.1 of Cohen's */
/*paper. Checked and it works, BUT with \p 50, it gives 17 correct digits   */
/*and \p 100 gives 26...                                                    */

\\ Define the parameters
p = 5;
k = 4

\\ Define the actual function
f(z) = (eta(z,1)*eta(5*z,1))^4; \\ f === delta5

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

R1(m) = {
    if(m <= p,
        I*intnum(t=1,[oo,2*Pi],V(I*t)*(m*I*t+1)^-k*f(I*t/(m*I*t+1))),
        I*intnum(t=1,[oo,2*Pi],V(I*t)*(I*t)^-k*f(-1/(I*t)))
    )
};


R2(m) = {
    if(m <= p,
        intnum(t=0,1,V(t+rho)*(m*(t+rho)+1)^(-k)*f((t+rho)/(m*(t+rho)+1))),
        intnum(t=0,1,V(t+rho)*(t+rho)^(-k)*f(-1/(t+rho)))
    );
};

\\ Index of Gamma_0(p) in SL_2(Z)
r = p + 1;

print("Computing..");

a(m) = {
    R1m = R1(m);
    R2m = R2(m);
    sum(n = 0, k-2,(-1)^n*binomial(k-2,n)*R1m[k-2-n+1]*conj(R2m[n+1]))
}

Peter_inner = sum(m = 1, r, a(m))/(2*I)^(k-1)/r;

print("The answer is :",precision(Peter_inner,50));
