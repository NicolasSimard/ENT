qhlinit(K) =
{
    my(hK = K.clgp.no, w, reps = redrepshnf(K), eiseval=vector(3,n,vector(hK)), tmp);
    w = if(imag(K.roots[1])>0,K.roots[1],conj(K.roots[1])); \\ make sure w in H
    
    \\ Evaluate the Eisenstein series at CM points    
    for(i=1,hK,
        tmp = subst(K.zk*reps[i],variable(K),w); \\ tmp = [a,(-b+sqrt(D))/2]
        eiseval[1][i] = tmp[1]^-2*G2star(tmp[2]/tmp[1]);
        eiseval[2][i] = tmp[1]^-4*G(4,tmp[2]/tmp[1]);
        eiseval[3][i] = tmp[1]^-6*G(6,tmp[2]/tmp[1]);
    );
    
    return([qhcinit(K),reps,eiseval]);
}
{
    addhelp(qhlinit,"qhlinit(K): initialize data to evaluate the L-function of a Hecke character of an imaginary quadratic field. This is a vector of the form [qhcinit(K),reps,eiseval].");
}

/* Tested: qhlfun(qhlinit(bnfinit('x^2+23)),[[i],[2,0]],2) = values in Watkin's paper

and

prod(i=0,2,qhlfun(qhlinit(bnfinit('x^2+23)),[[i],[2,0]],2))=
Pi^2/5*qhlfun(qhlinit(bnfinit('x^2+23)),[[0],[6,0]],4)

as in Watkin's paper.*/
qhlfun(qhldata,qhc,m) =
{
    \\ for the moment, we only accept T of the form [2*ell,0]
    if(qhc[2][2] != 0 || qhc[2][1]%2 == 1, error("Wrong values for infinity type."));
    
    my(t = qhc[2][1], K=qhldata[1][1], hK=K.clgp.no, reps = qhldata[2]);
    
    if(t == 0,
        if(vector(#qhc[1],i,qhc[1][i]%K.clgp[2][i]) == vector(#qhc[1]),
            error("L-function of K for trivial character has a pole at s=1.")
        );
        if(m != 1, error("Not implemented yet at s != 1."));
        return(-2*Pi/sqrt(abs(K.disc))*sum(i=1,#reps,qhceval(qhldata[1],qhc,reps[i])*log(sqrt(imag(idatouhp(K,reps[i])))*abs(eta(idatouhp(K,reps[i]),1))^2)))
    );
    
    if(2*m-t < 2 || t-m < 0, error("Not implemented for those values of m."));
    
    my(S, mf);
    
    mf = if(2*m-t == 2, delkformal('G2s,t-m), delkformal(GktoG4G6(2*m-t),t-m));
    S = sum(i=1,hK,
        qhceval(qhldata[1],qhc,qhldata[2][i])*subst(subst(subst(mf,'G2s,qhldata[3][1][i]),'G4,qhldata[3][2][i]),'G6,qhldata[3][3][i]);
    );
    return(I^t*(2*Pi)^m/gamma(m)*sqrt(abs(K.disc))^(t-m)*S);
}
{
    addhelp(qhlfun,"qhlfun(qhldata,qhc,m): Evaluates the Hecke L-function
    attached to the Hecke character [c,T] (where c are the components and T is
    the infinity type) at the point m (for the moment, we need t/2+1 <= m <= t, or t = 0 and m = 1).
    qhldata is the data returned by qhlinit.");
}
