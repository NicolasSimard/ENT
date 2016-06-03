\r ../lfunc/qhlfun-eisender.gp

/*Note that the petersson norm is normalized by removing the volume factor.*/
pnorm(qhdata,qhcomp,ell) =
{
    my(K=qhdata[1], hK = K.clgp.no);
    sqrt(abs(K.disc))*2*hK*(2*ell)!/(4*Pi)^(2*ell+1)*qhlfun(qhdata,qhcomp,4*ell,2*ell+2);
}