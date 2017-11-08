grammat(dim, f) = {
    my(M = matrix(dim, dim));
    for(i = 1, dim,
        for(j = i, dim,
            M[i,j] = f(i,j)));
    for(i = 2, dim,
        for(j = 1, i - 1,
            M[i,j] = conj(M[j,i])));
    M;
}

datatopipdata(data) = {
    if(type(data) == "t_INT", \\ data is a discriminant
        pipinit(bnfinit('x^2-data)),
        if(#data == 4, \\ data is a pipdata
            data,
            pipinit(data); \\ data is a number field?
        );       
    );
}

normalpipE2ell(pipdata, ell, ida, idb) = {
    if(ell <= 0, error("ell has to be > 0. Recieved ",ell));

    my(K = pipdata[1]);
    pip(pipdata,ell,ida,idb)\
    /E(2,idatolat(K,idealinv(K,ida)))^ell\
    /conj(E(2,idatolat(K,idealinv(K,idb))))^ell;
}

normalpipdnE2(pipdata, ell, ida, idb) = {
    if(ell <= 0, error("ell has to be > 0. Recieved ",ell));

    my(K = pipdata[1]);
    pip(pipdata,ell,ida,idb)\
    /dnE2(pipdata, ell-1, [idatolat(K,idealinv(K,ida))])[1]\
    /conj(dnE2(pipdata, ell-1, [idatolat(K,idealinv(K,idb))])[1]);
}

gramdetE2ell(data, ell, flag) = {
    if(ell <= 0, error("ell has to be > 0. Recieved ",ell));
    
    my(pipdata = datatopipdata(data), reps = redrepshnf(pipdata[1]), pol, M);
    
    M = grammat(#reps, (i,j) - > normalpipE2ell(pipdata, ell, reps[i], reps[j]));
            
    if(!flag, return(matdet(M)));
    
    pol = algdep(matdet(M),1);
    -polcoeff(pol,0)/polcoeff(pol,1);
}

gramdetdnE2(data, ell) = {
    if(ell <= 0, error("ell has to be > 0. Recieved ",ell));
    
    my(pipdata = datatopipdata(data), reps = redrepshnf(pipdata[1]), pol, M);
    
    matdet(grammat(#reps, (i,j) - > normalpipdnE2(pipdata, ell, reps[i], reps[j])));
}

gramdetpsi(data, ell) = {
    my(pipdata = datatopipdata(data), qhcs = qhchars(pipdata[1],[2*ell,0]));
    prod(i = 1, #qhcs, pnorm(pipdata, qhcs[i]));
}

normgramdetpsi(data, ell, Om, flag) = {
    my(pipdata = datatopipdata(data), Om_K, a, pol);
    Om_K = if(Om, Om, CSperiod(pipdata[1].disc));
    a = gramdetpsi(pipdata, ell)/Om_K^(4 * ell * pipdata[1].clgp.no);
    
    if(!flag, return(a));
    
    pol = algdep(a,1);
    -polcoeff(pol,0) / polcoeff(pol,1);
}

gramdetNida(data, ell) = {
    if(ell <= 0, error("ell has to be > 0. Recieved ",ell));
    
    my(pipdata = datatopipdata(data), reps = redrepshnf(pipdata[1]));
    
    matdet(grammat(#reps, (i,j) - > pip(pipdata, ell, reps[i], reps[j])))\
    /prod(i = 1, #reps, idealnorm(pipdata[1],reps[i]))^(2*ell);
    /* == matdet(matrix(#reps, #reps, i, j, pip(pipdata, ell, 1, idealmul(K, idealinv(K, reps[i]), reps[j]))));*/
}

normgramdetNida(data, ell, Om, flag) = {
    my(pipdata = datatopipdata(data), Om_K, a, pol);
    Om_K = if(Om, Om, CSperiod(pipdata[1].disc));
    a = gramdetNida(pipdata, ell)/Om_K^(4 * ell * pipdata[1].clgp.no);
    
    if(!flag, return(a));
    
    pol = algdep(a,1);
    -polcoeff(pol,0) / polcoeff(pol,1);
}

squarepart(N, ell) =  {
    my(M, sqpart);
    M = factor(inv);
    sqpart = prod(i=1,#M[,1],M[i,1]^(sign(M[i,2])*(abs(M[i,2])\2)));
    [inv/sqpart^2,sqpart];
}

/* issquareinvdenom(pipdata,ell=1) = issquare(denominator(invariant(pipdata,ell)));

factorinv(pipdata,ell=1) = factor(invariant(pipdata,ell));

factorinvdenom(pipdata,ell=1) = factor(denominator(invariant(pipdata,ell)));

factorinvnum(pipdata,ell=1) = factor(numerator(invariant(pipdata,ell)));

issquareinv(pipdata,ell=1) = issquare(invariant(pipdata,ell));*/