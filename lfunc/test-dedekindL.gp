/* This script contains test code for the dedekindL.gp script. The first test
is of course to compare it with the built-in Pari Dedekind zeta function. The
second is to compare the value of the automatically computed residue at s=1
with the one given by the class number formula.

There is a problem with the sign of the residue, however... The residue of the
L-function seems to be -1 times the residue given by the class number formula.
*/


read("dedekindL.gp");

print("--------------------Test 1------------------------");
fpol = x^5+2;

[DedekindZeta,DedekindFullGamma,DedekindRes] = initDedekindL(fpol);

bnf     = bnfinit(fpol);
LPari(s) = zetak(zetakinit(bnf),s);

\\Compare DedekindZeta(2) and built-in LPari(2)
print("L(2)                = ",DedekindZeta(2));
print(" (or using pari)    = ",LPari(2));

\\ Compare residue at 1 and the one given by class number formula
print("Residue at s=1");
print(" (automatically determined) = ",DedekindRes[1]/DedekindFullGamma(1));
residue = 2^bnf.r1*(2*Pi)^bnf.r2*bnf.reg*bnf.clgp.no/bnf.tu[1]/sqrt(abs(bnf.disc));
print(" (class number formula)     = ",residue);


print("--------------------Test 2------------------------");
fpol2 = x^3-6*x+3;

[DedekindZeta2,DedekindFullGamma2,DedekindRes2] = initDedekindL(fpol2);

print("Testing a theorem of Siegel");
print("(rationnal value at negative odd integers for totally real fields):");
print("L(-3)               = ",DedekindZeta2(-3));
print("L(-5)               = ",DedekindZeta2(-5));
print("L(-7)               = ",DedekindZeta2(-7));