K = bnfinit(x^2+47);
p = 11;
idp = idealprimedec(K,p)[1];
idpinv = idealinv(K,p);
idlat = idatolat(K,1);
idplat = idatolat(K,idpinv);
idplat2 = idatolat(K,idealpow(K,idpinv,2));
u1 = delta(idlat)/delta(idplat);
u2 = delta(idplat2)/delta(idplat);
u = u1^p*u2;