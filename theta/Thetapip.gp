{
    addhelp(Thetapip,"The Thetapip package is there to compute the Petersson
    inner product and norm of various theta functions attached to an imaginary
    quadratic field K. The first type of theta series is attached to an ideal
    of K and an integer ell>=0. The second is the theta series attached to a
    Hecke character of infinity type [2*ell,0], for some positive ell.
    
    Recall that a Hecke character of K is represented by [c,T].
    
    *Petersson inner product
    - pip(pipdata,ell,ida,idb) -> <theta_ida,theta_idb>
    - pnorm(data,qhc) -> <theta_qhc,theta_qhc>
    
    *Various objects attached to K
    - pipgrammat(pipdata,ell,{reps=redreps}) -> Gramm matrix in basis reps
    - pipgramdet(pipdata,ell,{reps=redreps}) -> determinant of pipgrammat
    - psigrammat(data,ell) -> Gramm matrix in basis of theta_psi
    - psigramdet(data,ell) -> determinant of psigrammat
    - transmat(K,ell,reps) -> transition matrix between reps and theta_psi
    
    *Auxilary functions
    - pipdatatoqhldata(pipdata) -> qhldata
    ");
}


/*-----------------Functions to compute to Petersson inner product -----------*/
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
    eiseval = matrix(hK,3,i,j,E(2*j,idatolat(K,reps[i])));
    
    return([K,reps,amb,eiseval]);
}
addhelp(pipinit,"pipinit(K,{verbose=0}): Initialise the data to compute the Petersson inner product of theta series attached to the imaginary quadratic field K. If verbose=1, display information about the quadratic field while the computations are performed.");

pnorm(data,qhc) =
{
    if(qhc[2][2] != 0, error("Wrong infinity type: ",qhc[2]));
    
    my(qhldata = if(#data == 4,pipdatatoqhldata(data),data));
    my(K=qhldata[1][1], hK = K.clgp.no, t = qhc[2][1]);
    
    sqrt(abs(K.disc))*2*hK*(t)!/(4*Pi)^(t+1)*qhlfun(qhldata,2*qhc,t+1);
}
addhelp(pnorm,"pnorm(data,qhc): return the norm of theta_psi, where psi is determined by qhc = [c,[2*ell,0]] and data is either the data returned by qhlinit or pipinit.");

pip(pipdata,ell,ida,idb) =
{
    if(ell == 0, error("ell has to be > 0."));
    
    my(K=pipdata[1], idbbar, coords, c0coords, c0, amb = pipdata[3], lambdac);
    idbbar = idealmul(K,idealinv(K,idb),idealnorm(K,idb));
    
    \\ Check if ida*idbbar is a square and return 0 if it isn't
    coords = bnfisprincipal(K,idealmul(K,ida,idbbar),0);
    for(i=1,#K.clgp.cyc,if(coords[i]%2 != 0 && K.clgp.cyc[i]%2 == 0, return(0)));
    
    \\ find an ideal c0 st ida*idbbar*c0^2 = lambdac0 O_K
    c0coords = vector(#coords,i,if(coords[i]%2==0,coords[i]/2,(coords[i]-K.clgp.cyc[i])/2));
    c0 = idealinv(K,idealfactorback(K,K.gen,c0coords));
    lambdac0 = generator(K,idealmul(K,ida,idealmul(K,idbbar,idealpow(K,c0,2))));
    
    \\ Compute the sum
    4*(abs(K.disc)/4)^ell*sum(i=1,#amb,(lambdac0*amb[i][2])^(2*ell)*d2l_1E2(pipdata,ell,idealmul(K,c0,amb[i][1])));
}
addhelp(pip,"pip(pipdata,ell,ida,idb,{flag=0}): Return the Petersson inner product of the theta series attached to ida and idb, with parameter ell. pipdata is the data returned by pipinit. By default, the Petersson inner product is normalized by removing the volume factor. If flag = 1, return the product without the constant C_K.");

\\ Evaluates d^(2*ell-1)E_2 at quadratic ideals
d2l_1E2(pipdata,ell,ida) =
{
    my(mu, K=pipdata[1], i0=1, coords = bnfisprincipal(pipdata[1],ida,0));
    \\i0 = select((rep->rep[2] == coords), pipdata[2],1)[1];
    while(pipdata[2][i0][2] != coords, i0 += 1);
    mu = generator(K,idealmul(K,idealinv(K,pipdata[2][i0][1]),ida)); \\ ida = mu*pipdata[2][i0][1]
    mu^(-4*ell)*substvec(delkformal('E2,2*ell-1),['E2,'E4,'E6],pipdata[4][i0,]);
}

/*-----------------------Various objects attached to K------------------------*/
psigrammat(data,ell) =
{
    my(K = if(#data == 4, data[1], data[1][1]), chars = qhchars(K,[2*ell,0]));
    matrix(K.clgp.no,K.clgp.no,i,j,if(i==j,pnorm(data,chars[i])));
}
addhelp(psigrammat,"psigrammat(data,ell): Gramm matrix of the Petersson inner product in the basis of theta_psi and data is either returned by qhlinit of pipinit.");

psigramdet(data,ell) = matdet(psigrammat(data,ell));

pipgrammat(pipdata,ell,reps) =
{
    my(hK = pipdata[1].clgp.no);
    if(reps == 0, reps = vector(hK,i,pipdata[2][i][1]));
    matrix(hK,hK,i,j,pip(pipdata,ell,reps[i],reps[j]));
}
addhelp(pipgrammat,"pipgrammat(pipdata,ell,{reps=redreps}): Gramm matrix of the Petersson inner product in the basis of reps (reduced reps by default.)");

pipgramdet(pipdata,ell,reps) = matdet(pipgrammat(pipdata,ell,reps));

transmat(K,ell,reps) =
{
    my(qhcdata = qhcinit(K), hK = K.clgp.no, ClK, chars = qhchars(K,[2*ell,0]));
    if(reps == 0, ClK = redrepshnf(K), ClK = parirepshnf(K));
    matrix(hK,hK,i,j,2/hK*qhceval(qhcdata,chars[i],ClK[j]));
}
addhelp(transmat,"transmat(K,ell,{reps=redreps}): Transition matrix M between the basis theta_psi and reps basis. It is such that M~*psigrammat*conj(M) = pipgrammat.");

/*-------------------------Auxilarry functions--------------------------------*/
pipdatatoqhldata(pipdata) = qhlinit(pipdata[1]);

generator(K,ida) = {
    my(tmp=bnfisprincipal(K,ida));
    if(tmp[1] != vector(#tmp[1])~, error("Non-principal ideal."));
    subst(K.zk*tmp[2],variable(K),K.roots[1]);
}
