/* Computing the symmetric square L-function of theta series.*/

/*Setting up the field.*/
ell = 2;
D = -47;
K = bnfinit(x^2-D);
hK = K.clgp.no;
wK = 2;

/*Setting up the modular form*/
k = 2*ell+1;
N = abs(D);
chi_D(n) = kronecker(D,n);

comps = [9];
qhc = [comps,[2*ell,0]];
anf = (n -> Vec(bintheta(K,qhc,1,'q,n)));

/*Setting up the symmetric square L-function*/
a(n) = {
    my(an=anf(n), ap2);
    direuler(p=2,n,
        if(N%p == 0,
            1/(1-p^(k-1)*X) \\ Not as given by Shimura's formula, but works
        ,
            ap2 = an[p]^2-chi_D(p)*p^(k-1);
            1/(1 - chi_D(p)*ap2*X + chi_D(p)*ap2*p^(k-1)*X^2 - p^(3*(k-1))*X^3)
        )
    )
};

astar = 1;
Vga = [0,1,2-k];
w = 2*k-1;
cond = N^2;
eps = 1;
L = lfuncreate([a,astar,Vga,w,cond,eps]);

printf("Estimated error on symmetric square L-function: %.2f%%", (1.0-abs(lfuncheckfeq(L))/default(realbitprecision))*100);

/* To have a nice functional equation, one must add Euler factors at the bad
primes to the Symmetric L-function defined in the thesis. At s=k, those Euler
factor comtribute to a factor of
\prod_{p|N}(1-p^-1)^-1 = N/eulerphi(N)
to the L-function. Hence the cancellation in the formula
Petersson norm = (Pi/2*eulerphi(N)*(4*Pi)^k/N^2/(k-1)!)^-1*lfun(L,k)
*/
(k-1)!*N/(4*Pi)^k*2/Pi*lfun(L,k)