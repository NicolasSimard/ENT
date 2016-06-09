\r ../Quadratic.gp
\r ../Modform.gp
\r Qhc.gp

qhlinit(K) =
{
    my(hK = K.clgp.no, w, reps, eiseval=vector(3,n,vector(hK)), tmp);
    w = if(imag(K.roots[1])>0,K.roots[1],conj(K.roots[1])); \\ make sure w in H
    
    reps = vector(hK,i,qfbtohnf(reduced_forms(K.disc)[i])); \\ Use idealred instead
    
    \\ Evaluate the Eisenstein series at CM points    
    for(i=1,hK,
        tmp = subst(K.zk*reps[i],variable(K),w); \\ tmp = [a,(-b+sqrt(D))/2]
        eiseval[1][i] = tmp[1]^-2*G2star(tmp[2]/tmp[1]);
        eiseval[2][i] = tmp[1]^-4*G(4,tmp[2]/tmp[1]);
        eiseval[3][i] = tmp[1]^-6*G(6,tmp[2]/tmp[1]);
    );
    
    return([qhcinit(K),reps,eiseval]);
}

/* Tested: qhlfun(qhcinit(-23),[i],2,2) = values in Watkin's paper
and prod(i=0,2,qhlfun(qhcinit(-23),[i],2,2))=Pi^2/5*qhlfun(qhcinit(-23),[0],6,4)
as in Watkin's paper.*/
qhlfun(qhldata,qhc,m) =
{
    \\ for the moment, we only accept T of the form [2*ell,0]
    if(qhc[2][2] != 0 || qhc[2][1]%2 == 1 || qhc[2][1] < 2 || 2*m-qhc[2][1] < 2 || qhc[2][1]-m < 0, error("Wrong values for infinity type and m."));
    
    my(t = qhc[2][1], S, mf, K=qhldata[1][1], hK=K.clgp.no);
    
    mf = if(2*m-t == 2, delkformal('G2s,t-m), delkformal(GktoG4G6(2*m-t),t-m));
    S = sum(i=1,hK,
        qhceval(qhldata[1],qhc,qhldata[2][i])*subst(subst(subst(mf,'G2s,qhldata[3][1][i]),'G4,qhldata[3][2][i]),'G6,qhldata[3][3][i]);
    );
    return(I^t*(2*Pi)^m/gamma(m)*sqrt(abs(K.disc))^(t-m)*S);
}
{
    addhelp(qhlfun,"qhlfun(qhdata,qhcomp,t,m): Evaluates the Hecke L-function
    attached to the Hecke character with components qhcomp and infinity type t
    at the point m (for the moment, we need t/2+1 <= m <= t). qhdata is the
    data returned by qhcinit.");
}