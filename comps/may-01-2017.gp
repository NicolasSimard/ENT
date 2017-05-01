/*Trying to see if kappa_psi is algebraic, since Baker's theorem doesn't apply*/

discs = discofclassno(5);

K = bnfinit(x^2+47);
reps = redrepshnf(K);
print(">>>K has discriminant ",K.disc);

hchar = [[1],[0,0]];
qhcdata = qhcinit(K);
psichar_1(ida) = qhceval(qhcdata,hchar,ida);

print(">>>Computing sum of psi^2: ");
sum(i=1,#reps,-psichar_1(reps[i])^2)

print(">>>Computing kappa psi: ");
kappa_psi=prod(i=1,#reps,(phi(K,reps[i])/phi(K,reps[1]))^(-psichar_1(reps[i])^2));

print(">>>Runing algdep on it (saved as f):");
print(f=algdep(kappa_psi,10));

/*print(">>>Reducing it (saved as fred):");
print(fred = polredbest(f));

print(">>>Calling quadhilbert:");
print(polredbest(quadhilbert(K.disc)));*/