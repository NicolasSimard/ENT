/*Try to attach a quantity to a pair of ideals ida and idb that only depends on
their ideal class.*/

K = bnfinit(x^2+47);
ell = 1;
T=[2*ell,0];
reps = redrepshnf(K);

p3=idealprimedec(K,3)[1];
p11=idealprimedec(K,11)[1];
ida=idealmul(K,p3,p11);

pipdata = pipinit(K);
N(comp,ell=ell) = pnorm(pipdata,[comp,[2*ell,0]]);
\\OmK = CSperiod(K.disc);

P(ida,idb,ell)=pip(pipdata,ell,ida,idb);

normalpip(ida,idb,ell)=
{
    P(ida,idb,ell)/E2star(idatolat(K,idealinv(K,ida)))^ell/conj(E2star(idatolat(K,idealinv(K,idb))))^ell;
}

normalpip2(ida,idb,ell)=
{   \\ This is simply P(...)/E2*((ida*\bar{idb})^-1)^(2*ell), but it is NOT an invariant!
    P(ida,idb,ell)/G(4,idatolat(K,idealmul(K,idealmul(K,idealinv(K,ida),idb),idealnorm(K,idb)^-1)))^ell;
}

inv1(ell)=matdet(matrix(#reps,#reps,i,j,normalpip(reps[i],reps[j],ell)));

inv2(ell)=matdet(matrix(#reps,#reps,i,j,pip(pipdata,ell,reps[i],reps[j])))/prod(i=1,#reps,idealnorm(K,reps[i]))^(2*ell);

\\ NOT an invariant...
inv3(ell,reps=reps) = matdet(matrix(#reps,#reps,i,j,normalpip2(reps[i],reps[j],ell)));

bernprimes(N) = {
    L=Set([]);
    for(n=2,N\2,L=Set(concat(L,factor(bernfrac(2*n))[,1]~)));
    L;
}

factorinv1(ell) = my(tmp=Vec(algdep(inv1(ell),1))); factor(-tmp[2]/tmp[1]);
