D = -7;
qhc = [[],[2,0]];
qhc2 = 2*qhc;
ell = qhc[2][1]/2;

K = bnfinit(x^2+7);
Ldata = qhcLdata(K,qhc2);

printf("Estimated error on symmetric square L-function: %.2f", lfuncheckfeq(lfuncreate(Ldata)));

lfun(lfuncreate(Ldata),2*ell+1)*lfun(lfuncreate(D),1)

4*K.clgp.no/K.tu[1]*sqrt(abs(D))*(2*ell)!/(4*Pi)^(2*ell+1)*lfun(lfuncreate(Ldata),2*ell+1)*3/8/Pi