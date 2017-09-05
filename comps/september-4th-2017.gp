/*An attempt to compute the Shimura-Maass operator on E2 differently.*/

Pr(r,z) = -(r==0)/24+suminf(n=1,n^r*sigma(n)*exp(2*Pi*I*n*z));

nhpart(r,z) = r!/(8*Pi*imag(z)^(r+1))/(4*Pi)^r;

pochhammer(a,m) = prod(i = 0, m-1, a + i);

c(n,z) = n!*(-1)^n/2/(4*Pi*imag(z))^(n+1);

s(f,n,z) = sum(r=0, n, (-1)^(n-r)*binomial(n,r)*pochhammer(2+r,n-r)/(4*Pi*imag(z))^(n-r)*f(r,z),0.);

f(n,z) = substvec(delkformal('E2,n),['E2,'E4,'E6],[E(2,z),E(4,z),E(6,z)]);	

f1(n,z) = c(n,z) +s(Pr,n,z); \\ == F(n,z)

f2(n,z) = s(nhpart,n,z) +s(Pr,n,z); \\ == F(n,z)

f3(n,z) = {
    c(n,z) 
    - (-1)^n*pochhammer(2,n)/(4*Pi*imag(z))^n/24
    + suminf(m=1,sum(r=0,n,(-1)^(n-r)*binomial(n,r)*pochhammer(2+r,n-r)/(4*Pi*imag(z))^(n-r)*m^r,0.)*sigma(m)*exp(2*Pi*I*m *z));
} \\ == F(n,z)

f4(n,z) = {
    (-1)^n*(1/(8*Pi*imag(z))-(n+1)/24)*n!/(4*Pi*imag(z))^n
    + suminf(m=1,sum(r=0,n,(-1)^(n-r)*binomial(n,r)*pochhammer(2+r,n-r)/(4*Pi*imag(z))^(n-r)*m^r,0.)*sigma(m)*exp(2*Pi*I*m *z));
} \\ == F(n,z)


F(K,n,ida) = my(L=idatolat(K,ida)); L[1]^-(2*n+2)*f(n,L[2]/L[1]); 

F1(K,n,ida) = my(L=idatolat(K,ida)); L[1]^-(2*n+2)*f1(n,L[2]/L[1]); 

F2(K,n,ida) = my(L=idatolat(K,ida)); L[1]^-(2*n+2)*f2(n,L[2]/L[1]); 

F3(K,n,ida) = my(L=idatolat(K,ida)); L[1]^-(2*n+2)*f3(n,L[2]/L[1]);

F4(K,n,ida) = my(L=idatolat(K,ida)); L[1]^-(2*n+2)*f4(n,L[2]/L[1]);
