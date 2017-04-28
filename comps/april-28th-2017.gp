/*Computing invariants of IQFs.*/

D = -23;
K = bnfinit('x^2-D);
pipdata = pipinit(K);

for (ell=1,15, print (ell,"::",factorinvdenom(pipdata,ell)));