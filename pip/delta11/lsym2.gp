/* Compute the Petersson norm of delta_11 using a special value of the
Symmetric square L-function attached to it. Note that delta5 is a newform.*/

k = 2; \\ Weight of delta11
N = 11;
anf(n) = lfunan(lfunetaquo([1,2;11,2]),n);

avec(n) = {
    my(an=anf(n));
    direuler(p=2,n,
        if(N%p==0,
            1/(1-an[p]^2*X), \\factor at the bad primes
            1/(1-(an[p]^2-p^(k-1))*X+p^(k-1)*(an[p]^2-p^(k-1))*X^2-p^(3*(k-1))*X^3)
        )
    )
};

a = (n -> avec(n));
astar = 1;
Vga = [0,1,2-k];
w = 2*k-1;
cond = N^2;
eps = 1;

Ldata = [a,astar,Vga,w,cond,eps];

printf("Estimated error on symmetric square L-function: %.2f%%", (1.0-abs(lfuncheckfeq(lfuncreate(Ldata)))/default(realbitprecision))*100);

print("The Petersson norm of delta11 is (stored as delta11_norm): ");
delta11_norm = (Pi/2*(4*Pi)^k/(k-1)!/N)^-1*lfun(lfuncreate(Ldata),k)/gammanot(N)
