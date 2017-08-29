/*Comparing the different algorithms to compute the Petersson norm of theta series*/

K = bnfinit(x^2+23);
ell = 2;
qhc = [[1],[2*ell,0]];

pipdata = pipinit(K);
Lsym2 = lfuncreate(lsym2data(K,qhc));
Lpsi2 = lfuncreate(qhcLdata(K,2*qhc));

print("With standard function: ",real(pnorm(pipdata,qhc)));
print("With L(psi^2,s):        ",real(4*K.clgp.no/K.tu[1]*sqrt(abs(K.disc))*(2*ell)!/(4*Pi)^(2*ell+1)*lfun(Lpsi2,2*ell+1)));
print("With Lsym2(s):          ",real((2*ell)!*abs(K.disc)/(4*Pi)^(2*ell+1)*2/Pi*lfun(Lsym2,2*ell+1)));
