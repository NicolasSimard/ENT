/*Compute the L-function attached to a Weber character modulo m of the number
field K, i.e. a character of the ray class group modulo m.

For example, the following cases were tested up to \p 50 and with different
component vectors. When the character was trivial, the result matched the
built-in zetak function (but it was considerably slower...):
K = bnfinit(y^2+23);
K = bnfinit(y^2+47);
K = bnfinit(y^2-1254);
K = bnfinit(y^3-175);

Input: fpol           -the polynomial defining the number field
       components     -the components of the character
       m              -modulus, i.e. an integral ideal and possibly som real places
*/

fpol       = y^3-199;
K          = bnfinit(fpol);
components = [2];

P13 = idealprimedec(K,13)[1];
\\m = [P13,[1,1]];  \\No archimedian places
m = [1,[1]];

\\ 1)--------------------Define the number field------------------------------
Kr = bnrinit(K,m,1);    \\ flag = 1 to have generators.

\\ 2)-------------------Compute the class group-------------------------------
Clm = Kr.clgp;           \\ [h_m,[d_1,...,d_n],[g_1,...,g_n]]

\\ 3)-----------------Define the character on Clm-----------------------------
/*A Weber character is defined by a component vector [a_1,...,a_n] of
#Clm.gen integers 0<=a_i<d_i. Then if the ideal A = g_1^e_1...g_n^e_n in
Clm, we define wChar(A) = exp(2*Pi*I*sum(i=1,n,a_i*e_i/d_i)).
*/

if(#components != #Clm.gen,error("The number of components is not good: ",\
    "Clm has ",#Clm.gen," generators (",#components," comp given)."));

wCharGen(charComp,A) = {
    local(e);
    e = bnrisprincipal(Kr,A,0); \\flag = 0 means we dont need to principalisation
    exp(2*Pi*I*sum(i=1,#Clm.gen,charComp[i]*e[i]/Clm.cyc[i]));
}

\\ Not sure but I think two ideals are coprime if their norm is.
isCoprime(K,A,B) = if(gcd(idealnorm(K,m[1]),idealnorm(K,A)) == 1,1,0);

wChar(A) = if(isCoprime(K,m,A), wCharGen(components,A), 0);

\\ 4)-------------Compute the coefficients of the Dirichlet series------------
/* The L-function of a Hilbert character (or Grossencharacter) is usually
defined via an Euler product. We use this Euler product to compute the
coefficients.*/

eulerFac(wChar,p) = {
    local(P);
    P = idealprimedec(K,p);
    prod(n=1,#P,1-wChar(P[n])*X^P[n].f);
}

a = direuler(p=2,10000,eulerFac(wChar,p)^(-1));

\\ 5)------------------Initialise the L-function------------------------------
print("\n\nInitialising the L-function of the number field of discriminant ",\
    K.disc," for the character of components ",components,"\n\n");

read("computel");

\\The exponential factor is sqrt(abs(K.disc)N(m_f)/Pi^[K:Q])^s
conductor = abs(K.disc)*idealnorm(K,m[1]);

\\Same as the L-function of a number field, but with a shift for the real places
nbInf = sum(i=1,#m[2],m[2][i]);
gammaV  = concat(vector(K.r1+K.r2-nbInf,x,0),vector(K.r2+nbInf,x,1));

\\Should be OK
weight = 1;

\\The sign will be determined if the character is not trivial. Otherwise 1.
\\sgn = if(sum(i=1,#components,components[i]),X,1);
sgn = X;

\\No poles as long as the character is not trivial
Lpoles    = if(sum(i=1,#components,components[i]),[],[1]);

initLdata("a[k]",,"conj(a[k])");

\\determine the sign
\\print(checkfeq());
sgneq = Vec(checkfeq());       \\ checkfeq() returns c1+X*c2, should be 0
\\sgn   = if(sum(i=1,#components,components[i]),-sgneq[2]/sgneq[1],1);
sgn = -sgneq[2]/sgneq[1];

\\ 6)---------------Do a few verifications------------------------------------
print("Sign          we will try to solve from the functional equation");
print("            = ", sgn);
print("|sign|      = ", abs(sgn));

print("Verifying functional equation. Error: ",errprint(checkfeq(1.1)));

\\ 7)-----------------Compute values of the L-function!-----------------------
/*\\ If the character is trivial, we compare our L-function with the built-in zetak
\\if(sum(i=1,#components,components[i]) == 0,\
\\print("zetak(2)    =         ",zetak(zetakinit(K),2)));

\\ If the character is non-trivial and K = Q(sqrt(-23)), we evaluate at  s = 1
\\ to compare with Watkin paper.
\\if(sum(i=1,#components,components[i]) && K.disc == -23,\
\\print("L(psi,1)    =         ",L(1)," as in watkins paper!"));
\\print("L(psi,2)    =         ",L(2));
\\print("  (check)   =         ",L(2,1.1));

\\ If the character is trivial, we compare the automatic evaluation of the
\\ residue at s = 1 with the formula obtained by the class number formula.
if(sum(i=1,#components,components[i]) == 0,{
print("Residue at s=1");
print(" (automatically determined) = ",Lresidues[1]);
\\ Determine the residue at s=1 using the class number formula
residue   = - 2^(K.r1+K.r2) * Pi^(K.r2/2) * K.reg * K.clgp.no / K.tu[1];
print(" (class number formula)     = ",residue);
});*/

print("testing values: ",L(2));
print("check           ",L(2,1.1));

