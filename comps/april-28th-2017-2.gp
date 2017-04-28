/*Testing Stark's example with all IQFs of class number 3...*/

minpol(K) = algdep(phi(K,redrepshnf(K)[2])/phi(K,1),3);

discs = discofclassno(3);

print(">>>Checking that the kappa_psi are units:");
for(i=1,#discs,print(discs[i],":",minpol(bnfinit('x^2-discs[i]))));

print(">>>Checking that the kappa_psi generate the Hilbert class field over K:");
for(i=1,#discs,print(discs[i],":\n",polredbest(minpol(bnfinit('x^2-discs[i]))),"\n",polredbest(quadhilbert(discs[i]))));
