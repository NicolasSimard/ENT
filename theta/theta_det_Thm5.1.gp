/*Computation of the theta determinant attached to a quadratic field.       */
/*The script theta_det_def.gp is too slow, so we use theorem 5.1 of Cohen's */
/*paper.                                                                    */

\\ Define the parameters
p = 23;
Darmonk = 1;
nbr_coeff = 50;

\\ Deffine the precision
\p 50;

print("Computing q-expansion");
\\ Compute the q-expansion and store it in a file as a list called qexps
output = concat("qexps/",concat(Str(p),concat("_",concat(Str(Darmonk),concat("_",concat(Str(nbr_coeff),".gp"))))));
commande = concat("python theta_qexp.py ",concat(Str(-p),concat(" ",concat(Str(Darmonk),concat(" ",concat(Str(nbr_coeff)," pari"))))));
system(concat(commande,concat(" > ",output)));

\\ Read the q-expansion. This defines qExps
read(output);

\\ The weight of the theta series
k = 2*Darmonk + 1;

\\ The class number of the quadratic field of discriminant -p
h = length(qExps);

\\ Define infinity and rho
oo = [1];
rho = exp(2*Pi*I/3);

f(i, z) = suminf(n=1, qExps[i][n]*exp(2*Pi*I*n*z));

V(z) = {
    local(v);
    v = [];
    for(i=0, k-2, v = concat(v,[z^i]));
    v
};

\\ We compute the integrals that go into the formula

R1(i,m) = {
    if(m > p,
        I*intnum(t=1,[oo,2*Pi],V(I*t)*(m*I*t+1)^-k*f(i,I*t/(m*I*t+1))),
        I*intnum(t=1,[oo,2*Pi],V(I*t)*(I*t)^-k*f(i,-1/(I*t)))
    )
};


R2(i,m) = {
    if(m > p,
        intnum(t=0,1,V(t+rho)*(m*(t+rho)+1)^(-k)*f(i,(t+rho)/(m*(t+rho)+1))),
        intnum(t=0,1,V(t+rho)*(t+rho)^(-k)*f(i,-1/(t+rho)))
    );
};

\\ Index of Gamma_0(p) in SL_2(Z)
r = p + 1;

Peter_inner(i,j) = {
    sum(m = 1, r,
        sum(n = 0, k-2,
            (-1)^n*binomial(k-2,n)*R1(i,m)[k-2-n+1]*conj(R2(j,m)[n+1])
        )
    )/(2*I)^(k-1)/r
};

M = matrix(h,h);

print("Computing the matrix...");
/*{for(i=1,h,
    for(j=1,i,
        M[i,j] = Peter_inner(i,j);
        M[j,i] = conj(M[i,j]);
        print((i*(i-1)/2+j)/(h*(h+1)/2)*100.,"%");
    );
);}*/

print(Peter_inner(2,2));
print("The answer is :",matdet(M));
