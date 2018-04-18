qhcs = qhchars(K,[2*ell,0]);
for(i = 1, K.clgp.no, print(algdep((pnorm(pipdata,qhcs[i])/Om_K^(4*ell))^K.clgp.no,2*K.clgp.no, round(0.9 * default(realprecision))),"\n"));