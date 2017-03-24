default(realprecision,200);

rankinprod(k,N,a,b) = {
    vector(min(#a,#b),n,
        sumdiv(n,m,if(issquare(n/m),2^omega(gcd(m,N))*real(a[m]*b[m])/m^(k-1)))
    )
};

w = 2;
N = 11;

conductor = N^2;
gammaV    = [0,1,w-1,w];
weight    = 1;
sgn       = 1;
Lpoles    = [1];

gammafactor(s) = N^s*Pi^(-2*s)*gamma(s/2)*gamma((s+1)/2)*gamma((s+w-1)/2)*gamma((s+w)/2);

/*
default(seriesprecision,cflength()+10);
an = Vec(bintheta(bnfinit('x^2+7),[[],[2,0]],1,'q));
bn = an;
*/

/*
an = vector(cflength()+5,n,ramanujantau(n));
bn = an;
*/

default(seriesprecision,cflength()+10);
an = Vec('q*eta('q)^2*eta('q^11)^2);
bn = an;

rankinprodcoeff = rankinprod(w,N,an,bn);
initLdata("rankinprodcoeff[k]");

print("Error in func. eq.  = ",errprint(checkfeq()));

tmp = 1; forprime(p=2,N,if(N%p == 0,tmp *= (p+1)/p));
pip = 6/Pi^2*gamma(w)/(4*Pi)^w/tmp*Lresidues[1]/gammafactor(1);
print(" Petersson norn (with Cohen's normalization) = ",pip/(3/Pi));