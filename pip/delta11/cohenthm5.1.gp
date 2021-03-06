/*Computation of the norm of delta5.                                        */
/*The script theta_det_def.gp is too slow, so we use theorem 5.1 of Cohen's */
/*paper. Checked and it works, BUT with \p 50, it gives 6 correct digits    */
/*and \p 100 gives 13 and \p 200 gives 23 and \p 400 gives 40.              */

\\ Define the parameters
p = 11;
k = 2;

\\ Define the actual function
f(z) = (eta(z,1)*eta(11*z,1))^2; \\ f === delta11

\\ Define rho
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
        I*intnum(t=1,[[1],2*Pi],V(I*t)*(m*I*t+1)^-k*f(I*t/(m*I*t+1))),
        I*intnum(t=1,[[1],2*Pi],V(I*t)*(I*t)^-k*f(-1/(I*t)))
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

Peter_norm = {
    sum(m = 1, r,
        sum(n = 0, k-2,
            (-1)^n*binomial(k-2,n)*R1(m)[k-2-n+1]*conj(R2(m)[n+1])
        )
    )/(2*I)^(k-1)/r
};

print("Computing..");
print("The answer is :",precision(Peter_norm,50));
