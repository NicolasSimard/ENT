{
    addhelp(Ellipticunits,"This package defines functions to compute units in the Hilbert class field H of imaginary quadratic fields K.
    
    Our main reference for the definitions is Iwasawa Theory for CM Elliptic
    curves of De Shalit.
    
    *Elliptic Units
    - fundtheta(L,z) -> z in C;
    - bigtheta(K,ellida,ida,z) -> z in C;
    
    *Siegel units
    - phi(K,ida) -> z in C;
    - deltaquot(K,ida,idb) -> z in C;
    - siegelunit(K,ida,idb) -> z in C;    
    ");
}

myeta(L,z) = Pi/area(L)*conj(z)+s2(L)*z;
addhelp(myeta,"myeta(L,z): Return the eta function attached to the lattices L in C evaluated at z: Pi/area(L)*conj(z)+s2(L)*z.");

fundtheta(L,z) = (2*Pi*I)^12*delta(L)*exp(-6*myeta(L,z)*z)*ellsigma(L,z)^12;
addhelp(fundtheta,"fundtheta(L,z): Return the theta function attached to the lattices L in C evaluated at z: (2*Pi*I)^12*delta(L)*exp(-6*myeta(L,z)*z)*ellsigma(L,z)^12.");

bigtheta(K,ellida,ida,z) = fundtheta(idatolat(K,ellida),z)^idealnorm(K,ida)/fundtheta(idatolat(K,idealmul(K,idealinv(K,ida),ellida)),z);
addhelp(bigtheta,"bigtheta(K,ellida,ida,z): Return the function big theta attached to fractional ideals ellida and ida of K evaluated at z. Big theta is just a twist of fundtheta(ellida,z) by N(ida)-sigma_ida. In general, ellida could be any lattice with CM by O_K, but to compute the twist, it is simpler to assume that this lattice is contained in K. This is the function used to generate elliptic units in ray class fields of K.");

phi(K,ida) = 
{
    my(L=idatolat(K,ida));
    sqrt(idealnorm(K,ida))*abs(L[1]^-1*eta(L[2]/L[1],1)^2);
}
addhelp(phi,"phi(K,ida): function N(ida)^1/2*|eta(ida)^2|.");


deltaquot(K,ida,idb=1) = delta(idatolat(K,idb))/delta(idatolat(K,idealmul(K,idealinv(K,ida),idb)));
addhelp(deltaquot,"deltaquot(K,ida,idb=1): return delta(idb)/delta(idb*ida^-1)");

siegelunit(K,ida,idb=1) = 
{
    my(idbinv, genidaidb, genidb);
    idbinv = idealinv(K,idb);
    genidb = subst(K.zk*bnfisprincipal(K,idealpow(K,idbinv,K.clgp.no))[2],variable(K),K.roots[1]);
    genidaidb = subst(K.zk*bnfisprincipal(K,idealpow(K,idealmul(K,ida,idbinv),K.clgp.no))[2],variable(K),K.roots[1]);
    abs(deltaquot(K,ida,idb)^K.clgp.no*(genidaidb/genidb)^12);
}
addhelp(siegelunit,"siegelunit(K,ida,idb=1): Return the Siegel unit attached to the ideal ida of K. This unit is known to land in the Hilbert class field of K. It is defined as (|Delta(idb)/Delta(idb*ida^-1)|N(ida)^6)^h_K");

/* Here are a few examples:
-----------------------------------Example 0-----------------------------------
(09:59) gp > \r Ellipticunits.gp
(09:59) gp > K=bnfinit(x^2+47);
(09:59) gp > \p 500
   realprecision = 500 significant digits
(09:59) gp > K=bnfinit(x^2+47);
(09:59) gp > p2=idealprimedec(K,2)[1];
(10:00) gp > p3=idealprimedec(K,3)[1];
(10:00) gp > p5=idealprimedec(K,5)[1];
(10:00) gp > ida=idealmul(K,p2,p3);
(10:01) gp > idb=idealmul(K,idealpow(K,p3,3),p5);
(10:01) gp > algdep(siegelunit(K,ida),10)
%23 = x^10 + 133547734060849615*x^9 + 4909059575945186345635433254174920*x^8 + 2934773274353554741585162128997411255*x^7 + 467669999439427716039105576882206218085*x^6 - 937365787675835092984768638688453323377*x^5 + 467669999439427716039105576882206218085*x^4 + 2934773274353554741585162128997411255*x^3 + 4909059575945186345635433254174920*x^2 + 133547734060849615*x + 1
(10:01) gp > algdep(siegelunit(K,ida,idb),10)
%24 = x^10 + 133547734060849615*x^9 + 4909059575945186345635433254174920*x^8 + 2934773274353554741585162128997411255*x^7 + 467669999439427716039105576882206218085*x^6 - 937365787675835092984768638688453323377*x^5 + 467669999439427716039105576882206218085*x^4 + 2934773274353554741585162128997411255*x^3 + 4909059575945186345635433254174920*x^2 + 133547734060849615*x + 1

As expected, the Siegel unit attached to ida and idb is a root of the same
polynomial as the one attached to ida only, since the later is a twist of the
former byt the Froebenius of idb^-1 (see deShalit).
-----------------------------------Example 1-----------------------------------
(16:11) gp > \p 1000;
(16:11) gp > \r Ellunits.gp ;
(16:12) gp > K=bnfinit(x^2+11);
(16:12) gp > idm=6;
(16:12) gp > ida=idealprimedec(K,5)[1];
(16:12) gp > v=quadgen(-11)/6.0;
(16:13) gp > factor(algdep(bigtheta(K,1,ida,v),20))
%9 =
[                                                                                                              x 4]

[x^6 + 55805258*x^5 + 7567746800224943*x^4 - 37653340800119188*x^3 + 11343835782131632111*x^2 + 6505864874*x + 1 1]

\\ We clearly see that bigtheta(K,1,ida,v) is a unit of degree 6 over Q!

(16:14) gp > bnrinit(K,idm).clgp.no
%10 = 6
(16:14) gp > bnrclassno(K,idm)
%11 = 6

\\ And here we clearly see that it generates the ray class field mod 6.

-----------------------------------Example 2-----------------------------------
As the following example suggests, those units do not always generate the ray class field.

(16:29) gp > \p 1000;
(16:30) gp > \r Ellunits.gp
(16:30) gp > K=bnfinit(x^2+23);
(16:30) gp > ida=idealprimedec(K,7)[1];
(16:30) gp > idm=6;
(16:30) gp > v=quadgen(-23)/6.0;
(16:30) gp > factor(algdep(bigtheta(K,1,ida,v)^-1,30))
%8 =
[x^3 - 22033926814609148578187136013720292189730965480141650607973202671269978*x^2 - 293244007989589654473076296484029147*x - 1 1]

(16:31) gp > bnrclassno(K,idm)
%9 = 6

\\ So bigtheta(K,1,ida,v)^-1 cannot generate the ray class field mod 6.
\\ In fact, it generates the Hilbert class field over K

(16:31) gp > polredbest(%8[1,1])
%10 = x^3 - x - 1

-------------------------------------Remark------------------------------------
Note that there is no loss of generality in taking ellida=O_K=1 (=L, using deShalit's
notation), since taking other lattices would give Galois conjugate units (prop
2.4 of "deShalit - Iwasawa theory of elliptic curves with complex multiplication").
However, the best trick to simplify computations is to take ellida=idm and ida
coprime to idm. Then one can simply take 1 as a primitive idm-torsion point.

-----------------------------------Example 3-----------------------------------
Using this remark, we are actually able to find a generator for the ray class
field of K for the modulus 6:

(10:10) gp > \p 10000
   realprecision = 10018 significant digits (10000 digits displayed)
(10:10) gp > K=bnfinit(x^2+23);
(10:10) gp > idm = 6;
(10:10) gp > p5 = idealprimedec(K,5)[1];
(10:10) gp > factor(algdep(bigtheta(K,idm,p5,1)^-1,10))
%15 =
[x 4]

[x^6 + ... - 3349541833612284587985321440972805964753946128592806*x + 1 1]

\\ So we found a generator of the ray class field mod 6.

-----------------------------------Example 4-----------------------------------
We illustrate the remark in another case.
(17:21) gp > \p 4000
   realprecision = 4007 significant digits (4000 digits displayed)
(17:21) gp > \r Ellunits.gp
(17:21) gp > K=bnfinit(x^2+47);
(17:21) gp > idm=6;
(17:21) gp > bnrclassno(K,idm)
%6 = 10
(17:21) gp > ida=idealprimedec(K,5)[1];
(17:22) gp > factor(algdep(bigtheta(K,idm,ida,1),10))
%8 = [x^10 + ... -3080580698[...]*x + 1 1]
(17:22) gp > bnrclassno(K,idm)
%24 = 10

\\ Which is a unit! Note that the height of this polynomial is huge and the
\\ computations must be done at very high precision.
*/