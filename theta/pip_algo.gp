/* This script uses he formula for the Petersson norm of theta series in terms
of linear combinations of derivatives of the Eisenstein series E2.*/

\r ../Modform.gp

pip(pipdata,ell,ida,idb,flag) = 
{
    my(K = pipdata[1],tmp);
    tmp=idealnorm(K,idb)^(2*ell)*S(pipdata,ell,idealmul(K,ida,idealinv(K,idb)));
    if(flag == 0, return(tmp)); \\ C_K*tmp = (.,.)
    if(flag == 1, return(4*(abs(K.disc)/4)^ell*tmp)); \\ V*K.clgp.no*(|K.disc|/4)^ell*tmptmp = (.,.)
}

S(pipdata,ell,ida) =
{
    my(mu, K=pipdata[1], amb = pipdata[3], cyc=K.clgp.cyc, coords, sqroot, c0, ac0);
    coords = bnfisprincipal(K,ida,0);
    for(i=1,#cyc,if(coords[i]%2 != 0 && cyc[i]%2 == 0, return(0)));
    sqroot = vector(#coords,i,if(coords[i]%2==0,coords[i]/2,(coords[i]-cyc[i])/2));
    c0 = idealinv(K,idealfactorback(K,K.gen,sqroot));
    ac0 = subst(K.zk*bnfisprincipal(K,idealmul(K,idealpow(K,c0,2),ida))[2],'x,K.roots[1]);
    ac0^(2*ell)*sum(i=1,#amb, amb[i][2]^(2*ell)*d2l_1E2(pipdata,ell,idealmul(K,c0,amb[i][1])));
}

\\ Evaluates d^(2*ell-1)E_2 at quadratic ideals
d2l_1E2(pipdata,ell,ida) =
{
    my(mu, K=pipdata[1], i0=1, coords = bnfisprincipal(pipdata[1],ida,0));
    while(pipdata[2][i0][2] != coords, i0 += 1);
    mu = bnfisprincipal(K,idealmul(K,idealinv(K,pipdata[2][i0][1]),ida))[2];
    mu = subst(K.zk*mu,'x,K.roots[1]); \\ ida = mu*pipdata[2][i0][1]
    mu^(-4*ell)*subst(subst(subst(delkformal('E2s,2*ell-1),'E2s,pipdata[4][1][i0]),'E4,pipdata[4][2][i0]),'E6,pipdata[4][3][i0]);
}

orient(ida) = if(imag(ida[2]/ida[1]) > 0, ida, vector(2,i,ida[3-i]));

\\ Petersson inner product init
pipinit(D,verbose) =
{
    if(!isfundamental(D), error(D," is not a fundamental discriminant."));
    my(tmp, K=bnfinit('x^2-D),w,hK=K.clgp.no);
    my(reps = vector(hK), amb = [], fs, eiseval = vector(3,n,vector(hK)));
    
    fs = reduced_forms(D);
    w = if(imag(K.roots[1])>0,K.roots[1],conj(K.roots[1])); \\ make sure w in H
    
    for(i=1,#fs,
        reps[i] = [qfbtohnf(fs[i]),bnfisprincipal(K,reps[i][1],0)];
        tmp = bnfisprincipal(K,idealpow(K,reps[i][1],2));
        if(tmp[1] == 0,
            amb = concat(amb,[[reps[i][1],subst(K.zk*tmp[2],'x,w)]]);
        );
    );
    if(verbose,
        print("Class group:                   ",K.clgp);
        print("Size of 2-torsion:             ",#amb);
        print("Size of ClK^2:                 ",K.clgp.no/#amb);
        print("Representatives:               ",reps);
    );
    
    \\ Evaluate the Eisenstein series at CM points
    if(verbose,print("Evaluating the Eisenstein series..."));
    for(i=1,hK,
        tmp = subst(K.zk*reps[i][1],'x,w); \\ tmp = [a,(-b+sqrt(D))/2]
        eiseval[1][i] = tmp[1]^-2*E2star(tmp[2]/tmp[1]);
        eiseval[2][i] = tmp[1]^-4*E(4,tmp[2]/tmp[1]);
        eiseval[3][i] = tmp[1]^-6*E(6,tmp[2]/tmp[1]);
    );
    
    return([K,reps,amb,eiseval]);
}

pipgrammat(pipdata,ell,reps='red) =
{
    my(hk = pipdata[1].clgp.no)
    if(reps = 'red,  reps = redrepshnf(pipdata[1]));
    if(reps = 'pari, reps = parirepshnf(pipdata[1]));
    matrix(hk,hk,i,j,pip(pipdata,ell,reps[i],reps[j],1));
}

pipgramdet(pipdata,ell,reps='red) = matdet(pipgrammat(pipdata,ell,reps));

/*
Mredreps(pipdata,ell) =
{
    my(hk=pipdata[1].clgp.no);
    matrix(hk,hk,i,j,pip(pipdata,ell,pipdata[2][i][1],pipdata[2][j][1],1));
}

MParireps(pipdata,ell) =
{
    my(K=pipdata[1], Clk=K.clgp, reps=[]);
    forvec(e=vector(#Clk.cyc,i,[0,Clk.cyc[i]-1]),
        reps=concat(reps,[idealfactorback(K,Clk.gen,e)]);
    );
    matrix(Clk.no,Clk.no,i,j,pip(pipdata,ell,reps[i],reps[j],1));
}*/

minpolZag(D,ell) = algdep(pipgramdet(pipinit(D),ell)/CSperiod(D)^(4*bnfclassno(D)*ell),5);