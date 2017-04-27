/*Testing Stark's example with IQF of class number 5...*/

discs = discofclassno(6*5);

K = bnfinit(x^2-discs[#discs]);
reps = redrepshnf(K);
print(">>>K has discriminant ",K.disc);

hchar = [[11],[0,0]];
qhcdata = qhcinit(K);
psichar_1(ida) = qhceval(qhcdata,hchar,ida);

print(">>>Computed using the definition (/h_K): ");
u=prod(i=1,#reps,phi(K,reps[i])^(-psichar_1(reps[i])^2));

print(">>>Runing algdep on it (saved as f):");
f=factor(algdep(u,K.clgp.no));
print(f=if(type(f) != "t_POL",f[#f[,1],1]));

print(">>>Reducing it (saved as fred):");
print(fred = polredbest(f));

print(">>>Calling quadhilbert:");
print(polredbest(quadhilbert(K.disc)));