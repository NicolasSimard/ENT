normalpip(pipdata,ell,ida,idb)=
{
    if(ell <= 0, error("ell has to be > 0. Recieved ",ell));

    my(K=pipdata[1]);
    pip(pipdata,ell,ida,idb)/G2star(idatolat(K,idealinv(K,ida)))^ell/conj(G2star(idatolat(K,idealinv(K,idb))))^ell;
}

invariant(data,ell) =
{
    if(default(realprecision)<500,warning("Precision is low (<500)..."));

    if(ell <= 0, error("ell has to be > 0. Recieved ",ell));
    
    my(pipdata);
    if(type(data) == "t_INT", \\ data is a discriminant
        pipdata = pipinit(bnfinit('x^2-data)),
        if(#data == 4, \\ data is a pipdata
            pipdata = data,
            pipdata = pipinit(data); \\ data is a number field?
        );       
    );
    
    my(reps=redrepshnf(pipdata[1]),pol);
    
    pol=algdep(matdet(matrix(#reps,#reps,i,j,normalpip(pipdata,ell,reps[i],reps[j]))),1);
    -polcoeff(pol,0)/polcoeff(pol,1);
}

squarepartinv(pipdata,ell) = 
{
    my(M,inv,sqpart);
    inv = invariant(pipdata,ell);
    M=factor(inv);
    sqpart = prod(i=1,#M[,1],M[i,1]^(sign(M[i,2])*(abs(M[i,2])\2)));
    [inv/sqpart^2,sqpart];
}

issquareinvdenom(pipdata,ell=1) = issquare(denominator(invariant(pipdata,ell)));

factorinv(pipdata,ell=1) = factor(invariant(pipdata,ell));

factorinvdenom(pipdata,ell=1) = factor(denominator(invariant(pipdata,ell)));

factorinvnum(pipdata,ell=1) = factor(numerator(invariant(pipdata,ell)));

issquareinv(pipdata,ell=1) = issquare(invariant(pipdata,ell));