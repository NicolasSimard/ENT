/* Computing the symmetric square L-function of Delta_7.*/

k = 3; \\ Weight of Delta_7
N = 7;
anf(n) = lfunan(lfunetaquo([1,3;7,3]),n);

chi(n) = kronecker(-7,n);

/*Setting up the symmetric square L-function*/
a(n) = {
    my(an=anf(n), ap2);
    direuler(p=2,n,
        if(N%p == 0,
            1/(1-p^(k-1)*X),\\ Not as given by Shimura's formula, but works
            ap2 = an[p]^2-chi(p)*p^(k-1);
            1/(1 - chi(p)*ap2*X + chi(p)*ap2*p^(k-1)*X^2 - p^(3*(k-1))*X^3)
        )
    )
};

astar = 1;
Vga = [0,1,2-k];
w = 2*k-1;
cond = N^2;
eps = 1;
Lsym2data = [a,astar,Vga,w,cond,eps];

printf("Estimated error on symmetric square L-function: %.2f", lfuncheckfeq(lfuncreate(Lsym2data)));

lfun(lfuncreate(Lsym2data),k)

(Pi/2*eulerphi(N)*(4*Pi)^k/N^2/(k-1)!)^-1*(lfun(lfuncreate(Lsym2data),k)*6/7)

/*The factor 6/7 = (1-1/7) = (1-7^(3-1)*p^-s)|s=3 is the factor at the bad
primes that is missing, in the defintion of L(Sym^2 f,...) given in the thesis,
to have a nice functionnal equation.*/