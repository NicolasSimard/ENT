/*L-function of Dirichlet characters.*/

/* Define generalized bernoulli numbers to verify computations. */
bernchi(n,dchar) = {
	my(F = dchar[1].mod, chi = lfunan(lfuncreate(dchar),F),P = bernpol(n));
	F^(n-1)*sum(a=1,F,chi[a]*subst(P,'x,a/F));
}

/* Gauss sums to define root numbers. */
gausssum(dchar) = {
	my(F = dchar[1].mod, chi = lfunan(lfuncreate(dchar),F));
	sum(a=1,F,chi[a]*exp(2*Pi*I*a/F));
}

/* Initial parameters. */
f = 23;
G = idealstar(,f);                         \\ (Z/fZ)^*
\\chi = vector(#G.cyc,i,random(G.cyc[i])); \\Some random characters, hoping it is primitive
chi = [1];
dchar = [G,chi];
n = 11 ;

/* Compute with built-in lfun. */
print("Evbaluatign at 1-",n,":");
print("Built-in value:              ",lfun(dchar,1-n));  						 \\Value at 1-n with built-in L-function creator
print("Using Bernoulli numbers:     ",-bernchi(n,dchar)/n);

/*Initialize by hand. */
a = n -> lfunan(lfuncreate(dchar),n);    \\ Could also do it by hand
astar = 1;								 \\ Sentinel value to indicate as = conj(a)
v = ideallog(,-1,G); sgn = round(prod(i=1,#chi,exp(2*Pi*I*v[i]*chi[i]/G.cyc[i]))); \\ dchar(-1) = +-1
Vga = [(sgn==-1)];					 \\ [0,0] if chi is even (i.e. sgn == 1), [0,1] otherwise
k = 1;									 \\ weight: s <-> k-sgn
N = G.mod;							 \\ if dchar is primitive, its conductor will be f...
eps = gausssum(dchar)/I^(sgn==-1)/sqrt(G.mod);	\\ root number
Ldata = [a,astar,Vga,k,N,eps];		 \\ 6 component vector, since no poles

print("Initialized by hand:         ", lfun(lfuncreate(Ldata),1-n));
print("(Check functionnal equation):", lfuncheckfeq(lfuncreate(Ldata)));
print("\n\n");

/* Find root number numerically. */
incLdata = [a,astar,Vga,k,N,'X];
print("Automatic root number computed:    ", lfunrootres(lfuncreate(incLdata))[3]);
print("Exact root number using Gauss sums:", eps);