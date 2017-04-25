K = bnfinit(x^2+47);
ell = 1;
T=[2*ell,0];
reps = redrepshnf(K);

p3=idealprimedec(K,3)[1];
p11=idealprimedec(K,11)[1];
ida=idealmul(K,p3,p11);

pipdata = pipinit(K);
N(comp,ell=ell) = pnorm(pipdata,[comp,[2*ell,0]]);
OmK = CSperiod(K.disc);

P(ida,idb,ell)=pip(pipdata,ell,ida,idb);

F(ida,idb,ell)=P(ida,idb,ell)/E2star(idatolat(K,idealinv(K,ida)))^ell/conj(E2star(idatolat(K,idealinv(K,idb))))^ell;

M(ell)=matrix(#reps,#reps,i,j,F(reps[i],reps[j],ell));

bernprimes(N) = {
    L=Set([]);
    for(n=2,N\2,L=Set(concat(L,factor(bernfrac(2*n))[,1]~)));
    L;
}