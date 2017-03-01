/*L-function of Dirichlet characters.*/

/* Define generalized bernoulli numbers to verify computations. */
bernchi(n,dchar) = {
	my(F = dchar[1].mod, chi = lfunan(lfuncreate(dchar),F),P = bernpol(n));
	F^(n-1)*sum(a=1,F,chi[a]*subst(P,'x,a/F));
}

gausssum(dchar) = {
	my(F = dchar[1].mod, chi = lfunan(lfuncreate(dchar),F));
	sum(a=1,F,chi[a]*exp(2*Pi*I*a/F));
}

N = 23;
G = idealstar(,N);
\\chi = vector(#G.cyc,i,random(G.cyc[i])); \\Some random characters, hoping it is primitive
chi = [1];
dchar = [G,chi];

n = 11 ;
print("Built-in value:         ",lfun(dchar,1-n));  						 \\Value at 1-n with built-in L-function creator
print("Using Bernoulli numbers:",-bernchi(n,dchar)/n);

/*Try to initialize by hand.*/
a = n -> lfunan(lfuncreate(dchar),n);    \\ Could also do it by hand
astar = 1;								 \\ Sentinel value to indicate astar = conj(a)
v = ideallog(,-1,G); sgn = round(prod(i=1,#chi,exp(2*Pi*I*v[i]*chi[i]/G.cyc[i]))); \\ dchar(-1) = +-1
Vga = [(sgn==-1)];					 \\ [0,0] if chi is even (i.e. sgn == 1), [0,1] otherwise
k = 1;									 \\ weight: s <-> k-sgn
cond = G.mod;							 \\ if dchar is primitive, its conductor will be N...
w = gausssum(dchar)/I^(sgn==-1)/sqrt(G.mod);	\\ root number
Ldata = [a,astar,Vga,k,cond,w];		 \\ 6 component vector, since no poles

print("Initialized by hand: ", lfun(lfuncreate(Ldata),1-n));
print("Check functionnal equation: ",lfuncheckfeq(lfuncreate(Ldata)));

/* Find root number numerically */
w=lfunrootres(lfuncreate([a,astar,Vga,k,cond,'X]))[3];
incLdata = [a,astar,Vga,k,cond,'X];
print("Automatic root number computed:    ", lfunrootres(lfuncreate(incLdata))[3]);
print("Exact root number using Gauss sums:", w);