/* Computing the symmetric square L-function of Delta_7.*/

k = 3; \\ Weight of Delta_7
N = 7;
anf(n) = lfunan(lfunetaquo([1,3;7,3]),n);

chi(n) = kronecker(-7,n);

a(n) = {
    my(an=anf(n), apnorm);
    direuler(p=2,n,
        1/(1-p^(k-1)*'X)/(1 - chi(p)*(an[p]^2-2*chi(p)*p^(k-1))*'X + chi(p)^2*p^(2*(k-1))*'X^2)
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

6/Pi^2*gamma(k)/(4*Pi)^k*7/8*lfun(lfuncreate(Lsym2data),k)
