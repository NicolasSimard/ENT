K = bnfinit(x^2 + 18377);
qhcdata = qhcinit(K);
qhc = [[15,1,0],[8,0]];

L = ideallist(K,default(seriesprecision));

F1(qhc)='q*Ser(apply((L->sum(i=1,#L,qhceval(qhcdata,qhc,L[i]))),L),'q);
F2(qhc)=bintheta(K,qhc);

print("Precision: ",default(realprecision),". Error < ", vecmax(abs(Vec(F1(qhc)-F2(qhc)))))