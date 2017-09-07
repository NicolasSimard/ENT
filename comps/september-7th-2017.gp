dnE2p(K,n,p,ida) = {
    my(reps, mu, idx, lat, z, tmp, coefv, v);
    
    reps = redrepshnf(K);
    [mu, idx] = idarep(K, reps, ida);
    
    lat = idatolat(K, reps[idx]);
    z = lat[2] / lat[1];
    coefv = vector(n + 1, r ,(-1)^(n - r + 1) * binomial(n, r - 1) * prod(i = 0, n - r, 1 + r + i)/(4 * Pi * imag(z))^(n - r + 1));
    (mu * lat[1])^-(2 * n + 2) * suminf(m = 1, if(m%p != 0, coefv * vector(n + 1, r, m^(r - 1))~ * sigma(m) * exp(2 * Pi * I * m * z)));
}

K = bnfinit('x^2+11);
p = 3;
ida = qfbtohnf(Heegner_forms(K.disc,p^2)[1]);
ell0 = p - 1;
ell(n) = ell0 + p^n*(p-1)/2;
printf("The base point is ell0 = %i", ell0);

[Om_C, M] = canperiod(K,ida,1);

N(ell) = dnE2p(K, 2*ell - 1, p, ida)*(Om_C * 2 * Pi * I)^(4 * ell);

\\ at \p 10000, factor(algdep(N(ell(1))-N(ell(0)),18)) with ell0 = 1 gives a polynomial of degree 6...