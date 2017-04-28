/*Testing Stark's example with IQF of class number 5...*/

minpolkappapsi(K,comp) =
{
    my(u,reps,qhcdata,fact);
    qhcdata = qhcinit(K);
    reps = redrepshnf(K);
    u=prod(i=1,#reps,phi(K,reps[i])^(-qhceval(qhcdata,[comp,[0,0]],reps[i])^2));
    fact = factor(algdep(u,K.clgp.no));
}

isinhilbertclassfield(K,poly) =
{
    my(f);
    f=polredbest(polcompositum(quadhilbert(K.disc),K.pol)[1]);
    f==polredbest(polcompositum(f,poly)[1]);
}

discs = discofclassno(6);

print(">>>Checking that the kappa_psi are units:");
for(i=1,#discs,print(discs[i],":",minpolkappapsi(bnfinit('x^2-discs[i]),[1])));

print(">>>Checking that the kappa_psi are in H:");
for(i=1,#discs,my(K=bnfinit('x^2-discs[i]));print(discs[i],":",isinhilbertclassfield(K,minpolkappapsi(K,[1]))));
