linroot(pol) = -polcoeff(pol, 0)/ polcoeff(pol, 1);

prodF(D, F) = {
    my(K = bnfinit(x^2-D), reps = parirepshnf(K));
    prod(i = 1, #reps, F(K, reps[i]));
}

mfCMprod(D, f, k) = {
    my(K = bnfinit(x^2-D));
    prodF(D, (K, ida) -> idealnorm(K, ida)^(k/2) * abs(f(idatolat(K, ida))));
}

mfCMprodnorm(D, f, k, deg = 1) = {
    my(tmp);
    /*tmp = mfCMprod(D, f, k) /(CSperiod(D)*sqrt(2))^(k * qfbclassno(D));*/
    tmp = mfCMprod(D, f, k) * (12 * sqrt(abs(D)))^(k/2*qfbclassno(D))/(CSperiod(D)*sqrt(2))^(k * qfbclassno(D));
    /*mfCMprod(D, x -> delta(x), 12) /(CSperiod(D)*sqrt(2))^(12 * qfbclassno(D)) == 1*/;
    
    if(!deg,
        return(tmp),
        if(deg == 1, linroot(algdep(tmp, deg)), algdep(tmp, deg)));
}

/*for(k = 1, 20, print(k,":", mfCMprodnorm(D, x -> E(2*k, x), 2*k, 2)))*/

normP(D, flag = 1) = {
    my(tmp, pol);
    tmp = P(D) * (12 * sqrt(abs(D)))^(qfbclassno(D))/CSperiod(D)^(2 * qfbclassno(D));
    
    if(!flag, return(tmp));
    
    pol = algdep(tmp,1);
    -polcoeff(pol,0)/polcoeff(pol,1);
}