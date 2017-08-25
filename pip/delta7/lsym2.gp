/* Computing the symmetric square L-function of Delta_7.*/

k = 3; \\ Weight of Delta_7
N = 7;
anf(n) = lfunan(lfunetaquo([1,3;7,3]),n);

chi(n) = kronecker(-7,n);

/*Setting up the symmetric square L-function*/
avec(n) = {
    my(an=anf(n));
    direuler(p=2,n,
        if(N%p == 0,
            1/(1-chi(p)*p^(k-1)*X),
            ap2 = an[p]^2-chi(p)*p^(k-1);
            1/(1 - ap2*X + chi(p)*p^(k-1)*ap2*X^2 - chi(p)*p^(3*(k-1))*X^3)
        )
    )
};

a = (n -> avec(n));
astar = 1;
Vga = [0,1,2-k];
w = 2*k-1;
cond = N^2;
eps = 1;

Lsym2data = [a,astar,Vga,w,cond,1,0];
tmp = lfunrootres(lfuncreate(Lsym2data))
Lsym2data[6] = tmp[3];
Lsym2data[7] = [[3,tmp[1]*x^-1+O(x^0)]];
L = lfuncreate(Lsym2data)
Lraw = lfuncreate([a,astar,Vga,w,cond,1,[[k,'r*'x^-1+O('x^0)]]]);

printf("Estimated error on symmetric square L-function: %.2f", lfuncheckfeq(L));

(Pi/2*eulerphi(N)/N^2*(4*Pi)^k/(k-1)!)^-1*lfun(lfuncreate(-7),1)*tmp[1]*(3/Pi/8)
lfun(Lraw,k)
(Pi/2*eulerphi(N)/N^2*(4*Pi)^k/(k-1)!)^-1*lfun(lfuncreate(-7),1)*lfun(Lraw,k)*(3/Pi/8)