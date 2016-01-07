/*This script is there to test the dirichletL.gp script.

The first test is o verifiy that the defaull Dirichlet L-function (with no
parameters) is the Riemann zeta function by comparing its values with the
built-in PARI zeta function.

The second test is to define the L-function of a random character.

The third test is to look at the special values of the L-function attached to
the kronecker symbol of a quadratic field.

The last test is to verify that the Dedekind zeta function of a quadratic
field factors as a product of the Riemann zeta function and the L-function of
the corresponding kronecker symbol.
*/

read("dirichletL.gp");

print("------------------------Test 1----------------------------");

[zetaQ,zetaQFullGamma] = initDirichletL();

print("zeta(2)    = ",zetaQ(2));
print("  (check)  = ",zeta(2));
print("zeta(I)    = ",zetaQ(I));
print("  (check)  = ",zeta(I));

print("------------------------Test 2----------------------------");
P = 23;
val = 11;
prim = znprimroot(P);
chi(p) = if(p%P,exp(2*Pi*I*val*znlog(p,prim)/(P-1)),0);

[Lchi,LchiFullGamma] = initDirichletL(chi,P);

print("L(chi,2)     = ",Lchi(2));
print("(check)      = ",Lchi(2,1.1));
print("L(chi,I)     = ",Lchi(I));
print("(check)      = ",Lchi(I,1.1));

print("------------------Test 3 (need prec > 28)-----------------");
fpol = x^2-53;

K = bnfinit(fpol);
chi(n) = kronecker(K.disc,n);
[Lchi,LchiFullGamma] = initDirichletL(chi,abs(K.disc));

print("L(chi,2)     = ",Lchi(2));
print("(check)      = ",Lchi(2,1.1));
print("Value at 0   = ",Lchi(0));
{
if(chi(-1) == 1,
print("Has to be    = ",0);,
print("Has to be    = ",-1/abs(K.disc)*sum(r=1,abs(K.disc)-1,chi(r)*r));
);
}
print("Value at 1   = ",Lchi(1));
print("CNF          = ",2^K.r1*(2*Pi)^K.r2*K.clgp.no*K.reg/K.tu[1]/sqrt(abs(K.disc)));
{
if(K.disc > 0,
    print("Value at -1  = ",Lchi(-1));
    somme = 0;
    k = 0;
    while(k^2 < K.disc,
        if(k%2 == K.disc%2, somme += 2*sigma((K.disc-k^2)/4,1));
        k += 1;
    );
    print("Has to be    = ",-somme/5.0);
);
}

print("------------------------Test 4----------------------------");
fpol = x^2-23;

K = bnfinit(fpol);
chi(n) = kronecker(K.disc,n);

read("DedekindL.gp");
[Lchi,LchiFullGamma] = initDirichletL(chi,abs(K.disc));
[zetaK,zetaKFullGamma,zetaKRes] = initDedekindL(fpol);

ss = 5;
print("Lchi(2)*zetaQ(2) = ",Lchi(ss)*zeta(ss));
print("zetaK(2)         = ",zetak(zetakinit(K),ss));

