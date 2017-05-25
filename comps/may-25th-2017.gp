K = bnfinit(x^2+163);
H = bnfinit(polredbest(polcompositum(quadhilbert(K.disc),K.pol)[1]));
ida = 1;
M = idatolat(K,ida);
j0 = ellj(M[2]/M[1]);
E_M = if(j0 == 1728,ellinit([1,0]),ellinit([1,0,0,-36/(j0-1728),-1/(j0-1728)]));
M_H = ellperiods(E_M); \\ Lattice of E_M/H
Om = M_H[2]/M[1];

z_M_H = (2*Pi*I)^4*E(4,M_H); \\ E(k,q) is not defined over Q, but (2*Pi*I)^k*E(k,q) is
z_M = (2*Pi*I)^4*E(4,M)/Om^4;