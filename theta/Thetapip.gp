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

pnorm(pipdata,qhc,algo) =
{
    my(K = pipdata[1], ell = qhc[2][1]/2);
    if(qhc[2][2] != 0, error("Wrong infinity type: ",qhc[2]));
    if(qhcistrivial(K,2*qhc), error("Genus character given! ",qhc));
    
    my(reps = pipdata[2], qhcdata = qhcinit(K));
    if(ell == 0,
        return(-K.clgp.no/3/K.tu[2]* sum(i=1,#reps,qhceval(qhcdata,2*qhc,reps[i])* log(idealnorm(K,reps[i])^6* abs(delta(idatolat(K,reps[i])))))));
    
    \\ ell > 0
    if(algo == "lpsi2",
        my(Lpsi2 = lfuncreate(qhcLdata(K,2*qhc)));
        return(4*K.clgp.no/K.tu[1]*sqrt(abs(K.disc))*(2*ell)!/(4*Pi)^(2*ell+1)*lfun(Lpsi2,2*ell+1)));
    
    \\ We use derivatives of E2
    my(d2l_1E2v = dnE2(pipdata,2*ell-1,reps,algo));
    4*K.clgp.no/K.tu[1]^2*(abs(K.disc)/4)^ell * sum(i=1,#reps,qhceval(qhcdata,2*qhc,reps[i])*d2l_1E2v[i]);
}
addhelp(pnorm,"pnorm(data,qhc): return the norm of theta_psi, where psi is determined by qhc = [c,[2*ell,0]] and data the data returned by pipinit.");

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
    lambdac0 = complexgen(K,idealmul(K,ida,idealmul(K,idbbar,idealpow(K,c0,2))));
    
    \\ Compute the sum
    4*(abs(K.disc)/4)^ell*sum(i=1,#amb,(lambdac0*amb[i][2])^(2*ell)*d2l_1E2(pipdata,ell,idealmul(K,c0,amb[i][1])));
}
addhelp(pip,"pip(pipdata,ell,ida,idb,{flag=0}): Return the Petersson inner product of the theta series attached to ida and idb, with parameter ell. pipdata is the data returned by pipinit. By default, the Petersson inner product is normalized by removing the volume factor. If flag = 1, return the product without the constant C_K.");

dnE2(pipdata,n,idav,algo) = {
    my(K = pipdata[1]);
    if(algo == "qexp",
        my(latv,zv,tmpv);
        latv = vector(#idav,i,idatolat(K,idav[i]));
        zv = vector(#latv,i,latv[i][2]/latv[i][1]);
        tmpv = vector(#zv, i, \\Vectorize... no need to recompute the coeffs...
            (-1)^n*(1/(8*Pi*imag(zv[i]))-(n+1)/24)*n!/(4*Pi*imag(zv[i]))^n
            + suminf(m=1,sum(r=0,n,(-1)^(n-r)*binomial(n,r)*prod(i=0,n-r-1,2+r+i)/(4*Pi*imag(zv[i]))^(n-r)*m^r,0.)*sigma(m)*exp(2*Pi*I*m*zv[i])));
        return(vector(#latv,i,latv[i][1]^-(2*n+2)*tmpv[i])));
        
    if(algo == "qexpv",
        my(latv, zv, tmpv, coefv, v);
        latv = vector(#idav, i, idatolat(K, idav[i]));
        zv = vector(#latv, i, latv[i][2] / latv[i][1]);
        coefv = vector(n + 1, r ,(-1)^(n - r + 1) * binomial(n, r - 1) * prod(i = 0, n - r, 1 + r + i)/(4 * Pi)^(n - r + 1));
        tmpv = vector(#zv, i,
            v = vector(n + 1, r, coefv[r] / imag(zv[i])^(n - r + 1));
            (-1)^n*(1 / (8 * Pi * imag(zv[i])) - (n + 1) / 24) * n! / (4 * Pi * imag(zv[i]))^n
            + suminf(m = 1, v * vector(n + 1, r, m^(r - 1))~ * sigma(m) * exp(2 * Pi * I * m * zv[i])));
        return(vector(#latv, i, latv[i][1]^-(2 * n + 2) * tmpv[i])));
    
    \\ At this point, we now we will use the polynomials d^nE2 in C[E2,E4,E6]
    my(pol, mu, idx, reps);
    reps = vector(#pipdata[2], i, pipdata[2][i][1]);
    pol = if(type(algo) == "t_CLOSURE", algo(n), delkformal('E2, n));
    vector(#idav, j, 
        [mu, idx] = idarep(K, reps, idav[j]);
        mu^-(2 * n + 2) * substvec(pol, ['E2, 'E4, 'E6], pipdata[4][idx,]));
}
{
addhelp(dnE2, "dnE2(pipdata,n,ida,{algo}): compute d^nE_2(ida) using algorithm algo. There are options are:
    1) algo = 0 (default): Use diffop to compute d^2E2 exactly as a polynomial in E2,E4 and E6.
    2) algo = pol(n), type(pol) = \"t_CLOSURE\": call the function pol to obtain the polynomial d^nE2.
    3) algo = \"qexp\": Use the q-expansion of d^nE2, without vectorizing.
    4) algo = \"qexpv\": Use the q-expansion of d^nE2 with vectorization.");
}

\\ Evaluates d^(2*ell-1)E_2 at quadratic ideals
d2l_1E2(pipdata,ell,ida) =
{
    my(mu, K=pipdata[1], i0=1, coords = bnfisprincipal(pipdata[1],ida,0));
    \\i0 = select((rep->rep[2] == coords), pipdata[2],1)[1];
    while(pipdata[2][i0][2] != coords, i0 += 1);
    mu = complexgen(K,idealmul(K,idealinv(K,pipdata[2][i0][1]),ida)); \\ ida = mu*pipdata[2][i0][1]
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

