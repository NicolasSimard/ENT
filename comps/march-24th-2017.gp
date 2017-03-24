K = bnfinit(x^2+71); hK = K.clgp.no; wK=K.tu[1];
ell = 10;
comps = [9];
qhc = [comps,[2*ell,0]];
qhcbar = [-comps,[qhc[2][2],qhc[2][1]]];

L = qhcLdata(K,2*qhc,1);
Lbar = qhcLdata(K,2*qhcbar,1);
pipdata = pipinit(K,0);

print("Computed using L-functions: ",4*hK/wK*(2*ell)!*sqrt(-K.disc)/(4*Pi)^(2*ell+1)*lfun(L,2*ell+1));
print("Computed using our formula: ",pnorm(pipdata,qhc));