{
    addhelp(Thetafunc,"This package defines the binary theta functions attached to imaginary quadratic fields. More precisely, it defines functions to compute q-expansions.
    
    *compute q-expansions
    - bintheta(K,X,ell,var='q,prec=0,flag=0) -> power series    
    ");
}

bintheta(K,X,ell,var='q,prec,flag) = {
    if((type(X) == "t_VEC" && #X == 3) || type(X) == "t_QFI", return(binthetaqfb(X,ell,var,prec,flag)));
    if(type(X) == "t_VEC" && #X == 2, return(binthetaqhc(K,X,var,prec,flag)));
    return(binthetaida(K,X,ell,var,prec,flag));
}
{
    addhelp(bintheta,"bintheta(K,X,ell,{var='q},{flag=0}): Return the q-expansion of the binary theta series attached to K, X and ell with 'q replaced by var. Here, X can be a Hecke character of K (i.e. X = [[c_1,...,c_d],[2*ell,0]]), in which case ell is ignored and the function returns the q-expansion of the associated newform, or X can be an ideal of K or a binary quadratic form of discriminant K.disc, in which case the corresponding q-expansion is returned. If flag = 0, the coefficients are expressed exactly in terms of w = quadgen(K.disc). If flag != 0, the coefficients are complex, which makes the computations slightly faster.")
}

/*
binthetaqhc1(K,qhc,var='q,flag) = {
    my(L=ideallist(K,default(seriesprecision)), data = qhcinit(K));
    var*Ser(vector(#L,n,sum(N=1,#L[n],qhceval(data,qhc,L[n][N]))),var);
}*/

binthetaqhc(K,qhc,var='q,prec,flag) = {
    my(wK=K.tu[1],ClK=redrepshnf(K),data=qhcinit(K));
    sum(i=1,#ClK,qhceval(data,qhc,ClK[i])^-1*binthetaida(K,ClK[i],qhc[2][1]/2,var,prec,flag))/wK;
}
{
    addhelp(binthetaqhc,"binthetaqhc(K,qhc,var='q,prec,flag): Return the q-expansion of the binary theta series attached to K and qhc with 'q replaced by var. Here, qhc is a Hecke character of K (i.e. X = [[c_1,...,c_d],[2*ell,0]]). If flag = 0, the coefficients are expressed exactly in terms of w = quadgen(K.disc). If flag != 0, the coefficients are complex, which makes the computations slightly faster.");
}

binthetaida(K,ida,ell,var='q,prec,flag) = binthetaqfb(idatoqfb(K,ida),ell,var,prec,flag);
{
    addhelp(binthetaida,"binthetaida(K,ida,ell,var='q,prec,flag): Return the q-expansion of the binary theta series attached to K, ida and ell with 'q replaced by var. Here, ida is an ideal of K. If flag = 0, the coefficients are expressed exactly in terms of w = quadgen(K.disc). If flag != 0, the coefficients are complex, which makes the computations slightly faster.");
}

binthetaqfb(f,ell,var='q,prec,flag) = {
    my(D, M, n, Zbasis, coeffs, L);

    if(prec == 0, prec = default(seriesprecision));
    f = Vec(f); \\ make sure f is a vector (not of type t_QFI)
    D = f[2]^2-4*f[1]*f[3];
    coeffs = vector(prec+1);
    
    L = mateigen([f[1],f[2]/2;f[2]/2,f[3]],1)[1];
    M = floor(sqrt(prec/min(L[1],L[2])));
    
    Zbasis = [f[1],-(f[2] + D%2)/2+quadgen(D)];
    if(flag, Zbasis[2] *= 1.0);
    
    forvec(pt=[[-M,M],[-M,M]],
        n = f[1]*pt[1]^2 + f[2]*pt[1]*pt[2] + f[3]*pt[2]^2;
        if(n <= prec, coeffs[n+1] += (pt[1]*Zbasis[1]-pt[2]*Zbasis[2])^(2*ell));
    );
    
    subst(Ser(coeffs,'X),'X,var);
}
{
    addhelp(binthetaqfb,"binthetaqfb(f,ell,var='q,flag): Return the q-expansion of the binary theta series attached to f and ell with 'q replaced by var. Here, f is a positive definite quadratc form. If flag = 0, the coefficients are expressed exactly in terms of w = quadgen(K.disc). If flag != 0, the coefficients are complex, which makes the computations slightly faster.");
}
