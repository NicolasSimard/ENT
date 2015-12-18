/* This script uses the formula of Darmon-Lauder-Rotger to compute the
Petersson norm of the theta series attached to an imaginary quadratic
field for various Hecke characters.

We first do it for Hecke characters of infinity type [0,0], since these are
also Hilbert characters and so they were already implemented earlier.

Not sure this will work, since L(psi^2,s) has a pole at 1 when psi^2 = 1...
*/

P = 7;
disc = -P;               \\ Has to be of the form -p for the moment
k = 0;
if(disc%4 > 1, error("Not a discriminant!"));

K = bnfinit(x^2-disc);
ClK = K.clgp;            \\ [h_K,[d_1,...,d_n],[g_1,...,g_n]]


\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\\         Define the L-function of the Hecke characters                  \\\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

components = [];

if(#components != #ClK.gen,error("The number of components is not good: ",\
    "ClK has ",#ClK.gen," generators (",#components," comp given)."));

hCharGen(charComp,A) = {
    if(#charComp == 0, return(1));
    local(e);
    e = bnfisprincipal(K,A,0); \\flag = 0 means we dont need to principalisation
    exp(2*Pi*I*sum(i=1,#ClK.gen,charComp[i]*e[i]/ClK.cyc[i]));
}

hChar(A) = hCharGen(components,A);

eulerFac(hChar,p) = {
    local(P);
    P = idealprimedec(K,p);
    prod(n=1,#P,1-hChar(P[n])*X^P[n].f);
}

a = direuler(p=2,5000,eulerFac(hChar,p)^(-1));

print("\n\nInitialising the L-function of the number field of discriminant ",\
    K.disc," for the character of components ",components,"\n\n");

read("../computel");

conductor = abs(K.disc);  \\The exponential factor is sqrt(abs(K.disc)/Pi^2)^s
gammaV  = concat(vector(K.r1+K.r2,x,0),vector(K.r2,x,1)); \\Same as the L-function of a number field
weight = 1;               \\Should be OK
sgn = if(sum(i=1,#components,components[i]),X,1);
Lpoles = if(sum(i=1,#components,components[i]),[],[1]); \\No poles as long as the character is not trivial

initLdata("a[k]",,"conj(a[k])");

\\determine the sign
sgneq = Vec(checkfeq());       \\ checkfeq() returns c1+X*c2, should be 0
sgn   = if(sum(i=1,#components,components[i]),-sgneq[2]/sgneq[1],1);

peterNorm = 6*(2*k)!*K.clgp.no/(Pi*(4*Pi)^(2*k+1)*sqrt(abs(disc)))*\
            P/(P+1)*L(2*k+1);

print(L(2));
print(L(2,1.1));
