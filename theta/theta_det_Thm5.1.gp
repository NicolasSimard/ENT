/*Computation of the theta determinant attached to a quadratic field. The
script theta_det_def.gp is too slow, so we use theorem 5.1 of Cohen's paper.
*/

\rtheta_qexp.gp;

eval_form(coeff, z) = suminf(n=1, coeff[n]*exp(2*Pi*I*n*z));

V(z, k) = {
    local(v);
    v = [];
    for(i=0, k-2, v = concat(v,[z^i]));
    v
};

\\ We compute the integrals that go into the formula

int1(f,k,m,p) = {
    if(m <= p,
        I*intnum(t=1,[[1],2*Pi],V(I*t,k)*(m*I*t+1)^-k*eval_form(f,I*t/(m*I*t+1))),
        I*intnum(t=1,[[1],2*Pi],V(I*t,k)*(I*t)^-k*eval_form(f,-1/(I*t)))
    )
};


int2(f,k,m,p) = {
    local(rho);
    rho = exp(2*Pi*I/3);
    if(m <= p,
        intnum(t=0,1,V(t+rho,k)*(m*(t+rho)+1)^-k*eval_form(f,(t+rho)/(m*(t+rho)+1))),
        intnum(t=0,1,V(t+rho,k)*(t+rho)^-k*eval_form(f,-1/(t+rho)))
    );
};

Petersson_inner(f1, f2, k, p) = {
    local(Lim,Ljm,S);
    S = 0;
    for(m = 1, p + 1,
        int1m = int1(f1,k,m,p);
        int2m = int2(f2,k,m,p);
        S += sum(n = 0, k-2, (-1)^n*binomial(k-2,n)*int1m[k-2-n+1]*conj(int2m[n+1]))
    );
    S/(2*I)^(k-1)/(p + 1)
};

theta_det(p, Darmonk, verbose = 0) = {
    if(verbose, print("Computing q-expansion."));
    nbr_coeff = 100;
    qexps = theta_qexps(-p, Darmonk, nbr_coeff);

    \\ The weight of the theta series
    k = 2*Darmonk + 1;

    \\ The class number of the quadratic field of discriminant -p
    h = length(qexps);
    M = matrix(h,h);
    if(verbose, print("Computing the matrix (using thm 5.1)."));
    for(i=1,h,
        for(j=1,i,
            M[i,j] = Petersson_inner(qexps[i],qexps[j], k, p);
            M[j,i] = conj(M[i,j]);
            if(verbose, print(round((i*(i-1)/2+j)/(h*(h+1)/2)*100.),"%"));
        );
    );
    matdet(M)
};

