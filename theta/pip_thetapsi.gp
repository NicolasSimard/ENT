\r ../lfunc/qhlfun-eisender.gp

/*Note that the petersson norm is normalized by removing the volume factor.*/
pnorm(qhldata,qhc) =
{
    if(qhc[2][2] != 0, error("Wrong infinity type: ",qhc[2]));
    my(K=qhldata[1][1], hK = K.clgp.no, t = qhc[2][1], qhcsq);
    qhcsq = [vector(#qhc[1],i,2*qhc[1][i]),[2*t,0]];
    sqrt(abs(K.disc))*2*hK*(t)!/(4*Pi)^(t+1)*qhlfun(qhldata,qhcsq,t+1);
}
{
    addhelp(pnorm,"pnorm(qhdata,qhcomp,ell): return the norm of theta psi, where
    psi is determined by qhcomp and of infinity type 2*ell.");
}

psigrammat(qhldata,ell) =
{
    my(K = qhldata[1][1], chars = qhchars(K,[2*ell,0]));
    matrix(K.clgp.no,K.clgp.no,i,j,if(i==j,pnorm(qhldata,chars[i])));
}

psigramdet(qhldata,ell) = matdet(psigrammat(qhldata,ell));

transmat(qhcdata,ell,reps) =
{
    my(K = qhcdata[1], hK = K.clgp.no, ClK, chars = qhchars(K,[2*ell,0]));
    if(reps == 0, ClK = redrepshnf(K), ClK = parirepshnf(K));
    matrix(hK,hK,i,j,2/hK*qhceval(qhcdata,chars[i],ClK[j]));
}
{
    addhelp(transmat,"transmat(qhcdata,ell,reps): Transition matrix M between the
    basis theta_psi and theta_ida. It is such that
    M~*psigrammat*conj(M) = pipgrammat.");
}