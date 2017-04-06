/* Compute the Petersson norm of delta_11 using a special value of the
twisted (by the trivial character mod 11) Symmetric square L-function
attached to it, as in Hida's paper on congruences of cusp forms.
Note that delta11 is a newform. Doesn't seem to work...*/

k = 2; \\ Weight of delta11
N = 11;
anf(n) = lfunan(lfunetaquo([1,2;11,2]),n);

a_twist(n) = {
    my(an=anf(n));
    direuler(p=2,n,
        if(N%p==0,
            1, \\factor at the bad primes
            1/(1-(an[p]^2-p^(k-1))*X+p^(k-1)*(an[p]^2-p^(k-1))*X^2-p^(3*(k-1))*X^3)
        )
    )
};
astar = 1;
Vga = [0,1,2-k];
w = 2*k-1;
cond = N^2;
eps = 1;

Ldata = [a_twist,astar,Vga,w,cond,eps];

printf("Estimated error on symmetric square L-function: %.2f%%", (1.0-abs(lfuncheckfeq(lfuncreate(Ldata)))/default(realbitprecision))*100);

print("Using Hida's formula: ");
(Pi*(4*Pi)^k/(k-1)!)^-1*N*lfun(lfuncreate(Ldata),k)/gammanot(N)
