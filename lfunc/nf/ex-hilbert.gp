/*Compute the L-function attached to a non-trivial Hilbert character of a
number field K, i.e. a character of the class group.

For example, the following cases were tested up to \p 50 and with different
component vectors. When the character was trivial, the result matched the
built-in zetak function (but it was considerably slower...):
K = bnfinit(y^2+23); => K.clgp = [3,[3],...]
K = bnfinit(y^2+47); => K.clgp = [5,[5],...]
K = bnfinit(y^2-1254); => K.clgp = [4,[2,2],]
K = bnfinit(y^3-175); => K.clgp = [3,[3],...] \\Computes the root number well
K = bnfinit(y^3-19); => K.clgp = [3,[3],...] \\Computes the root number well
K = bnfinit(y^3-37); => K.clgp = [3,[3],...] \\Computes the root number well
K = bnfinit(y^3-199); => K.clgp = [9,[9],...] \\index 3

Note that a formula of Watkins says that the root number for Hilbert
characters is simply hChar(K.diff)^-1. In particular, if Z_K has a power
basis, then the different is principal (D_K = (f'(\alpha)) if f is the min
pol of \alpha) and so the root number is one. However, it seems that we must
not take the inverse in its formula...

The following code was used to find number fields of class number an index >1:
forprime(p=1,N,{K=bnfinit(y^3-p); print(p,": ",K.clgp.no,",",K.index)})

Note that by a theorem of Bardestani, the density of monogenic (i.e. index 1)
number fields in the familly K(y^q-p) is at least (q-1)/q wrt p. In
particular, at least 2/3 of the number fields above are monogenic.

Input: fpol           -the polynomial defining the number field
       components     -the components of the character
*/

fpol       = y^3-175;
components = [1];

\\ 1)--------------------Define the number field------------------------------
K = bnfinit(fpol);

\\ 2)-------------------Compute the class group-------------------------------
ClK = K.clgp;            \\ [h_K,[d_1,...,d_n],[g_1,...,g_n]]

\\ 3)-----------------Define the character on ClK-----------------------------
/* A Hilbert character is defined by a component vector [a_1,...,a_n] of
#ClK.gen integers 0<=a_i<d_i. Then if the ideal A = g_1^e_1...g_n^e_n in
ClK, we define hChar(A) = exp(2*Pi*I*sum(i=1,n,a_i*e_i/d_i)).
*/

if(#components != #ClK.gen,error("The number of components is not good: ",\
    "ClK has ",#ClK.gen," generators (",#components," comp given)."));

hCharGen(charComp,A) = {
    local(e);
    e = bnfisprincipal(K,A,0); \\flag = 0 means we dont need to principalisation
    exp(2*Pi*I*sum(i=1,#ClK.gen,charComp[i]*e[i]/ClK.cyc[i]));
}

hChar(A) = hCharGen(components,A);

\\ 4)-------------Compute the coefficients of the Dirichlet series------------
/* The L-function of a Hilbert character (or Grossencharacter) is usually
defined via an Euler product. We use this Euler product to compute the
coefficients.*/

eulerFac(hChar,p) = {
    local(P);
    P = idealprimedec(K,p);
    prod(n=1,#P,1-hChar(P[n])*X^P[n].f);
}

a = direuler(p=2,30000,eulerFac(hChar,p)^(-1));

\\ 5)------------------Initialise the L-function------------------------------
print("\n\nInitialising the L-function of the number field of discriminant ",\
    K.disc," for the character of components ",components,"\n\n");

read("computel");

\\The exponential factor is sqrt(abs(K.disc)/Pi^2)^s
conductor = abs(K.disc);

\\Same as the L-function of a number field
gammaV  = concat(vector(K.r1+K.r2,x,0),vector(K.r2,x,1));

\\Should be OK
weight = 1;

\\The sign will be determined if the character is not trivial. Otherwise 1.
sgn = if(sum(i=1,#components,components[i]),X,1);

\\No poles as long as the character is not trivial
Lpoles = if(sum(i=1,#components,components[i]),[],[1]);


initLdata("a[k]",,"conj(a[k])");

\\determine the sign
\\print(checkfeq());
sgneq = Vec(checkfeq());       \\ checkfeq() returns c1+X*c2, should be 0
sgn   = if(sum(i=1,#components,components[i]),-sgneq[2]/sgneq[1],1);

\\ 6)---------------Do a few verifications------------------------------------
print("Sign          we will try to solve from the functional equation");
print("                                 = ", sgn);
print("|sign|                           = ", abs(sgn));
print("Verify with a formula of Watkins = ", hChar(K.diff));

print("Verifying functional equation. Error: ",errprint(checkfeq(1.1)));

\\ 7)-----------------Compute values of the L-function!-----------------------
\\ If the character is trivial, we compare our L-function with the built-in zetak
if(sum(i=1,#components,components[i]) == 0,\
print("zetak(2)    =         ",zetak(zetakinit(K),2)));

\\ If the character is non-trivial and K = Q(sqrt(-23)), we evaluate at  s = 1
\\ to compare with Watkin paper.
if(sum(i=1,#components,components[i]) && K.disc == -23,\
print("L(psi,1)    =         ",L(1)," as in watkins paper!"));
print("L(psi,2)    =         ",L(2));
print("  (check)   =         ",L(2,1.1));

\\ If the character is trivial, we compare the automatic evaluation of the
\\ residue at s = 1 with the formula obtained by the class number formula.
if(sum(i=1,#components,components[i]) == 0,{
print("Residue at s=1");
print(" (automatically determined) = ",Lresidues[1]);
\\ Determine the residue at s=1 using the class number formula
residue   = - 2^(K.r1+K.r2) * Pi^(K.r2/2) * K.reg * K.clgp.no / K.tu[1];
print(" (class number formula)     = ",residue);
});

