/* Computing the symmetric square L-function of Delta_7.*/

k = 3; \\ Weight of Delta_7
N = 7;
anf(n) = lfunan(lfunetaquo([1,3;7,3]),n);

chi(n) = kronecker(-7,n);

/*Setting up the symmetric square L-function*/
a(n) = {
    my(an=anf(n));
    direuler(p=2,n,
		ap2 = an[p]^2-chi(p)*p^(k-1);
        1/(1 - chi(p)*ap2*X + chi(p)*ap2*p^(k-1)*X^2 - p^(3*(k-1))*X^3)
    )
};
astar = 1;
Vga = [0,1,2-k];
w = 2*k-1;
cond = N^2;
Lsym2data(cond,eps) = [a,astar,Vga,w,cond,eps];
{   
for(i=0,10,
    eps = lfunrootres(lfuncreate(Lsym2data(N^i,0)))[3];
    L = lfuncreate(Lsym2data(N^i,eps));
    print(lfuncheckfeq(L))
);
}
/*eps = lfunrootres(lfuncreate(Lsym2data))[3]
sym2data[6] = eps;
L = lfuncreate(Lsym2data);

printf("Estimated error on symmetric square L-function: %.2f", lfuncheckfeq(L));*/

(Pi/2*eulerphi(N)*(4*Pi)^k/N^2/(k-1)!)^-1*lfun(L,k)*(3/Pi/8)