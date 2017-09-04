/*An attempt to compute the Shimura-Maass operator on E2 differently.*/

Pr(r,z) = -(r==0)/24+suminf(n=1,n^r*sigma(n)*exp(2*Pi*I*n*z));

nhpart(r,z) = r!/(8*Pi*imag(z)^(r+1))/(4*Pi)^r;

pochhammer(a,m) = prod(i = 0, m-1, a + i);

c(n,z) = (-1)^n*sum(r=0,n,n!*binomial(n+1,r+1)/(r+1)/(2*I)^r/(4*Pi)^(n-r));

s(f,n,z) = sum(r=0, n, (-1)^(n-r)*binomial(n,r)*pochhammer(2+r,n-r)/(4*Pi*imag(z))^(n-r)*f(r,z));

del(f,n,z,k=2) = c(n,z)/(8*Pi*imag(z)^n+1) +s(f,n,z,k);

F0(n,z) = substvec(delkformal('E2,n),['E2,'E4,'E6],[E(2,z),E(4,z),E(6,z)]);	

F1(n,ida) = c(n,z)/(8*Pi*imag(z)^(n+1)) +s(Pr,n,z);

F2(n,z) = s(nhpart,n,z) +s(Pr,n,z); \\ == F0(ell,z)