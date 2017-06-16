K = bnfinit('x^2+23);
h_K = K.clgp.no;
ell = ellinit(ellfromj(ellj(idatouhp(K,1))));
L = ellperiods(ell);

print("\n\nFor E_2: ",factor(algdep((2*Pi*I)^2*E(2,L),2*h_K)));
print("\n\nFor E_4: ",factor(algdep((2*Pi*I)^4*E(4,L),2*h_K)));
print("\n\nFor E_6: ",factor(algdep((2*Pi*I)^6*E(6,L),2*h_K)));