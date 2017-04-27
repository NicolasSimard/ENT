/*Generalize Stark's example to all IQF of class number 3...*/

classno3 = discofclassno(3);

K = bnfinit(x^2-classno3[4]);
reps = redrepshnf(K);
C = reps[2]; \\ Representative of a non-trivial class
print(">>>K has discriminant ",K.disc);

hchar = [[1],[0,0]];
qhcdata = qhcinit(K);
psichar_1(ida) = qhceval(qhcdata,hchar,ida);

print(">>>Computed using the definition (/3): ");
print(prod(i=1,3,phi(K,reps[i])^(-psichar_1(reps[i])^2)));

print(">>>Computed using the shortcut (saved as u): ");
print(u=phi(K,C)/phi(K,1));

print(">>>Runing algdep on it (saved as f):");
print(f=algdep(u,3));

print(">>>Reducing it (saved as fred):");
print(fred = polredbest(f));

print(">>>Calling quadhilbert:");
print(polredbest(quadhilbert(K.disc)));