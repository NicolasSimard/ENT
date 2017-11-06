K = bnfinit(x^2+239);
reps = redrepshnf(K);

hchar = [[10],[0,0]];
qhcdata = qhcinit(K);
psichar(ida) = qhceval(qhcdata,hchar,ida);

\\print(vector(#reps,i,psichar(reps[i])^2));

print(">>>Computing kappa psi: ");
kappa_psi=prod(i=1,#reps,phi(K,reps[i])^(-psichar(reps[i])^2));

print(">>>Computing algdep of kappa_psi: ");
factor(algdep(kappa_psi,40))

