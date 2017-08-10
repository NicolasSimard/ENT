/* Computing the symmetric square L-function of Delta_7.*/

k = 3; \\ Weight of Delta_7
N = 7;
anf(n) = lfunan(lfunetaquo([1,3;7,3]),n);

chi(n) = kronecker(-7,n);

/*Setting up the symmetric square L-function*/
a(n) = {
    my(an=anf(n));
    direuler(p=2,n,
        if(N%p==0,
            1/(1-an[p]^2/p^(k-1)*X), \\factor at the bad primes
            1/(1-chi(p)*an[p]/p^((k-1)/2)*X+chi(p)*X^2)
        )
    )
};
astar = 1;
Vga = [1,k-1,k];
w = 1;
cond = N^2;
incLsym2data(cond) ={
    eps = lfunrootres(lfuncreate([a,astar,Vga,w,cond,0]))[3];
    [a,astar,Vga,w,cond,eps]
}

printf("Estimated error on symmetric square L-function: %.2f", lfuncheckfeq(lfuncreate(Lsym2data)));

lfun(lfuncreate(Lsym2data),1)

6/Pi^2*gamma(k)/(4*Pi)^k*7/8*lfun(lfuncreate(Lsym2data),1)