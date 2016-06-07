/* This script uses he formula for the Petersson norm of theta series in terms
of linear combinations of derivatives of the Eisenstein series E2.*/

\r ../Modform.gp
\r ../Quadratic.gp

pip(pipdata,ell,ida,idb,flag) = 
{
    my(K = pipdata[1],tmp);
    tmp=idealnorm(K,idb)^(2*ell)*S(pipdata,ell,idealmul(K,ida,idealinv(K,idb)));
    if(flag == 0, return(4*(abs(K.disc)/4)^ell*tmp)); \\ V^-1*4*(|K.disc|/4)^ell*tmp = (.,.)
    if(flag == 1, return(tmp)); \\ C_K*tmp = (.,.)
}

S(pipdata,ell,ida) =
{
    my(mu, K=pipdata[1], amb = pipdata[3], cyc=K.clgp.cyc, coords, sqroot, c0, ac0);
    coords = bnfisprincipal(K,ida,0);
    for(i=1,#cyc,if(coords[i]%2 != 0 && cyc[i]%2 == 0, return(0)));
    sqroot = vector(#coords,i,if(coords[i]%2==0,coords[i]/2,(coords[i]-cyc[i])/2));
    c0 = idealinv(K,idealfactorback(K,K.gen,sqroot));
    ac0 = subst(K.zk*bnfisprincipal(K,idealmul(K,idealpow(K,c0,2),ida))[2],variable(K),K.roots[1]);
    ac0^(2*ell)*sum(i=1,#amb, amb[i][2]^(2*ell)*d2l_1E2(pipdata,ell,idealmul(K,c0,amb[i][1])));
}

\\ Evaluates d^(2*ell-1)E_2 at quadratic ideals
d2l_1E2(pipdata,ell,ida) =
{
    my(mu, K=pipdata[1], i0=1, coords = bnfisprincipal(pipdata[1],ida,0));
    while(pipdata[2][i0][2] != coords, i0 += 1);
    mu = bnfisprincipal(K,idealmul(K,idealinv(K,pipdata[2][i0][1]),ida))[2];
    mu = subst(K.zk*mu,variable(K),K.roots[1]); \\ ida = mu*pipdata[2][i0][1]
    mu^(-4*ell)*subst(subst(subst(delkformal('E2s,2*ell-1),'E2s,pipdata[4][1][i0]),'E4,pipdata[4][2][i0]),'E6,pipdata[4][3][i0]);
}

orient(ida) = if(imag(ida[2]/ida[1]) > 0, ida, vector(2,i,ida[3-i]));

\\ Petersson inner product init
pipinit(K,verbose) =
{
    my(tmp,sq,w,hK=K.clgp.no);
    my(reps = [], amb = [], eiseval = vector(3,n,vector(hK)));
    
    w = if(imag(K.roots[1])>0,K.roots[1],conj(K.roots[1]));
    
    forvec(e=vector(#K.clgp.cyc,i,[0,K.clgp.cyc[i]-1]),
        tmp = idealred(K,idealfactorback(K,K.clgp.gen,e));
        reps = concat(reps,[[tmp,e~]]);
        sq = bnfisprincipal(K,idealpow(K,tmp,2));
        if(sq[1] == 0,
            amb = concat(amb,[[tmp,subst(K.zk*sq[2],variable(K),w)]]);
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
        tmp = subst(K.zk*reps[i][1],variable(K),w); \\ tmp = [a,(-b+sqrt(D))/2]
        eiseval[1][i] = tmp[1]^-2*E2star(tmp[2]/tmp[1]);
        eiseval[2][i] = tmp[1]^-4*E(4,tmp[2]/tmp[1]);
        eiseval[3][i] = tmp[1]^-6*E(6,tmp[2]/tmp[1]);
    );
    
    return([K,reps,amb,eiseval]);
}

pipgrammat(pipdata,ell,reps=0) =
{
    my(hK = pipdata[1].clgp.no, ClK);
    if(reps == 0, ClK = vector(hK,i,pipdata[2][i][1]));
    matrix(hK,hK,i,j,pip(pipdata,ell,ClK[i],ClK[j]));
}

pipgramdet(pipdata,ell,reps=0) = matdet(pipgrammat(pipdata,ell,reps));
