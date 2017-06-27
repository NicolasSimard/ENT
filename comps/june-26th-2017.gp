K = bnfinit('x^2+19);
p = 7;
ell = 3;
n = 2;
ell2 = ell+p^n*(p-1);

[Om_C, M] = canperiod(K,1,1);
pipdata = pipinit(K);

print("Evaluation at ",ell,":");
F_ell = pip(pipdata,ell,1,1)*(Om_C*2*Pi*I)^(4*ell);
M1 = factor(algdep(F_ell,2*K.clgp.no))

print("Evaluation at ",ell2,":");
F_ell2 = pip(pipdata,ell2,1,1)*(Om_C*2*Pi*I)^(4*ell2);
M2 = factor(algdep(F_ell2,2*K.clgp.no))

print("Difference: ");
M3 = factor(algdep(F_ell-F_ell2,2*K.clgp.no))