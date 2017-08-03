/* Computing the twisted by chi_D symmetric square L-function of theta series.*/

/*Setting up the field.*/
ell = 2;
D = -7;
K = bnfinit(x^2-D);
hK = K.clgp.no;
wK = 2;

/*Setting up the modular form*/
k = 2*ell+1;
N = abs(D);
chi_D(n) = kronecker(D,n);

comps = [];
qhc = [comps,[2*ell,0]];
anf = (n -> Vec(bintheta(K,qhc,1,'q,n)));

/*Setting up the symmetric square L-function*/
a(n) = {
    my(an=anf(n));
    direuler(p=2,n,
		if(N%p == 0,
		1,
		ap2 = an[p]^2-chi_D(p)*p^(k-1);
        1/(1-chi_D(p)*ap2*X+chi_D(p)*p^(k-1)*ap2*X^2-p^(3*(k-1))*X^3);
		)
    )
};
astar = 1;
Vga = [0,1,2-k];
w = 2*k-1;
cond = N^5;
eps = lfunrootres(lfuncreate([a,astar,Vga,w,cond,'X]))[3];
Lsym2data = [a,astar,Vga,w,cond,eps];

printf("Estimated error on symmetric square L-function: %.2f%%", (1.0-abs(lfuncheckfeq(lfuncreate(Lsym2data)))/default(realbitprecision))*100);
/*
\\defining L(qhc^2,s)
Lqhc2data = qhcLdata(K,2*qhc,1);

printf("Estimated error on Hecke L-function: %.2f%%", (1.0-abs(lfuncheckfeq(lfuncreate(Lqhc2data)))/default(realbitprecision))*100);

\\define zeta_N(s)
zeta_N(N,s) = {
	my(ps = factor(N)[1..-1,1]);
	prod(i=1,#ps,(1-ps[i]^-s))*zeta(s);
}

s=10;
print("Checking Lemma 4 of the notes:");
lfun(lfuncreate(Lsym2data),s)-lfun(lfuncreate(Lqhc2data),s)*zeta_N(N,s+1-k)

print("Computed using L-functions: ");
\\((4*Pi)^k/(k-1)!*eulerphi(N)/2/hK)^-1*lfunrootres(lfuncreate(Ldata))[1];*/