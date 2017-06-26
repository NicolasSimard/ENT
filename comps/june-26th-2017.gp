K = bnfinit('x^2+11);
p = 5;

ell = 3;
n = 1;
ell2 = ell+p^n*(p-1);

[Om_C, M] = canperiod(K,ida,1);
pipdata = pipinit(K);
print("Evaluation at ",ell,":");
F_ell = pip(pipdata,ell,1,1)*(Om_C*2*Pi*I)^(4*ell)
factor(algdep(F_ell,2*K.clgp.no))
f1 = Vec(factor(algdep(F_ell,2*K.clgp.no))[1,1]);
r1 = -f1[2]/f1[1];
print("Evaluation at ",ell2,":");
F_ell2 = pip(pipdata,ell2,1,1)*(Om_C*2*Pi*I)^(4*ell2)
factor(algdep(F_ell2,2*K.clgp.no))
f2 = Vec(factor(algdep(F_ell2,2*K.clgp.no))[1,1]);
r2 = -f2[2]/f2[1];

