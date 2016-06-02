\r ../Quadratic.gp
\r ../Modform.gp

/* Given a fundamental discriminant D, creates the imaginary quadratic field
K and computes the principalisations y_i of the generators of ClK, i.e. if g_i is
a generator of ClK of order o_i, then g_i^o_i=y_i O_K.*/
qhcinit(D) =
{
    if(!isfundamental(D) || D >= -4, error(D," is not a fundamental discriminant or D >= -4."));
    my(mu, K=bnfinit('x^2-D),ClK = K.clgp, princ=vector(#ClK.cyc));
    for(i=1,#ClK.cyc,
        mu = bnfisprincipal(K,idealpow(K,ClK.gen[i],ClK.cyc[i]))[2];
        princ[i] = subst(K.zk*mu,'x,K.roots[1]);
    );
    [K,princ];
}

/* Any Hecke character of K evaluated at a generator g_i of the class group has
value psi(g_i) = y_i^(t/o_i)*zeta_i^c_i, where y_i is as above, t is the
infinity type, o_i is the order of g_i, zeta_i = exp(2*Pi*I/o_i) and 0 <=c_i<o_i.
Then if ida = mu*g_1^e_1*...*g_d^e_d, one can determine psi(ida) from the above
information.*/
qhchar(qhdata,qhcomp,t,ida) =
{
    my(K=qhdata[1],ClK = K.clgp,decomp = bnfisprincipal(K,ida));
    if(#qhcomp != #ClK.cyc, error("Invalid component vector."));
    mu = subst(K.zk*decomp[2],'x,K.roots[1]);
    mu^t*prod(i=1,#ClK.cyc,(qhdata[2][i]^(t/ClK.cyc[i])*exp(2*Pi*I*qhcomp[i]/ClK.cyc[i]))^decomp[1][i]);
}

/* Tested:
qhlfun(qhcinit(-23),[i],2,2) = values in Watkin's paper.
*/
qhlfun(qhdata,qhcomp,t,m) =
{
    if(t%2 == 1 || t < 2 || 2*m-t < 2 || t-m < 0, error("Wrong values for infinity type and m."));
    
    my(S, mf, tmp, K=qhdata[1], w, hK=K.clgp.no, eiseval=vector(3,i,vector(hK)));
    
    fs = reduced_forms(K.disc);
    w = if(imag(K.roots[1])>0,K.roots[1],conj(K.roots[1])); \\ make sure w in H
    
    reps = vector(hK,i,qfbtohnf(fs[i]));
    
    \\ Evaluate the Eisenstein series at CM points    
    for(i=1,hK,
        tmp = subst(K.zk*reps[i],'x,w); \\ tmp = [a,(-b+sqrt(D))/2]
        eiseval[1][i] = tmp[1]^-2*E2star(tmp[2]/tmp[1]);
        eiseval[2][i] = tmp[1]^-4*E(4,tmp[2]/tmp[1]);
        eiseval[3][i] = tmp[1]^-6*E(6,tmp[2]/tmp[1]);
    );
    
    mf = if(2*m-t == 2, delkformal('E2s,t-m), delkformal(EktoE4E6(2*m-t),t-m));
    S = sum(i=1,hK,
        qhchar(qhdata,qhcomp,t,reps[i])*subst(subst(subst(mf,'E2s,eiseval[1][i]),'E4,eiseval[2][i]),'E6,eiseval[3][i]);
    );
    return(I^t*(2*Pi)^m/gamma(m)*sqrt(abs(K.disc))^(t-m)*S);
}