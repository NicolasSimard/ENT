kloosterman_sum(chi,N,m,n,c) = \
    sum(r=1,c,if(gcd(r,c) == 1,chi(r)*exp(2*Pi*I*(m*r+n*lift(Mod(1/r,c)))/c),0));

poincarre_series_qexp(m,chi,N,k,bound) = {
    local(coeffs,temp);
    coeffs = [];
    for(n=1,bound,
        \\ c = N*l otherwise the Kloosterman sum is 0
        temp = suminf(l=1,kloosterman_sum(chi,N,m,n,N*l)*besselj(k-1,4*Pi*sqrt(m*n)/(N*l))/(N*l));
        coeffs = concat(coeffs,[kronecker(n,m)^((k-1)/2)*2*Pi/I^k*temp]);
    );
    coeffs;
}
