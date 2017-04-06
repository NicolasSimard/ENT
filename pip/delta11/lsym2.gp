/* Compute the Petersson norm of delta_11 using a special value of the
Symmetric square L-function attached to it. Note that delta5 is a newform.*/

w = 2; \\ Weight of delta11
N = 11;

avec(n) = {
    my(an=lfunan(lfunetaquo([1,2;11,2]),n));
    direuler(p=2,n,
        if(N%p==0,
            1/(1-an[p]^2*X), \\factor at the bad primes
            1/(1-(an[p]^2-p^(w-1))*X+p^(w-1)*(an[p]^2-p^(w-1))*X^2-p^(3*(w-1))*X^3)
        )
    )
};

a = (n -> avec(n));
astar = 1;
Vga = [0,1,2-w];
k = 2*w-1;
cond = N^2;
eps = 1;

printf("Estimated error: %.2f%%", (1.0-abs(lfuncheckfeq(lfuncreate([a,astar,Vga,k,cond,eps])))/default(realbitprecision))*100);

/*
index = N;
forprime(p=2,N,if(N%p == 0,index = index*(1+p^-1)));

norme = ((4*Pi)^w/factorial(w-1)*Pi/N^2/eulerphi(N))^-1/index*L(w);
other_norm = 2/Pi*factorial(w-1)*N/(4*Pi)^w*index*L(1);
print("The norm should be   = ",norme);
print("or                   = ",other_norm);*/

