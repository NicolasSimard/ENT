\r ../lfunc/qhlfun-eisender.gp

/*Note that the petersson norm is normalized by removing the volume factor.*/
pnorm(qhdata,qhcomp,ell) =
{
    my(K=qhdata[1], hK = K.clgp.no);
    sqrt(abs(K.disc))*2*hK*(2*ell)!/(4*Pi)^(2*ell+1)*qhlfun(qhdata,qhcomp,4*ell,2*ell+1);
}
{
    addhelp(pnorm,"pnorm(qhdata,qhcomp,ell): return the norm of theta psi, where
    psi is determined by qhcomp and of infinity type 2*ell.");
}

transmat(qhdata,ell,reps) =
{
    my(K = qhdata[1], hK = K.clgp.no, ClK, qhcomps=[]);
    if(reps == 0, ClK = redrepshnf(K), ClK = reps);
    forvec(e=vector(#K.clgp.cyc,i,[0,K.clgp.cyc[i]-1]),
        qhcomps = concat(qhcomps,[e]);
    );
    matrix(hK,hK,i,j,2/hK*qhchar(qhdata,qhcomps[i],2*ell,ClK[j]));
}
{
    addhelp(transmat,"transmat(qhdata,ell,reps): Transition matrix M between the
    basis theta_psi and theta_ida. It is such that
    det((M^*)*psigrammat*M) = det(pipgrammat).");
}