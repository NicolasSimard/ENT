/*Analyse difference between CM elliptic curves defined by ideals that differ
by an ambiguous ideal, i.e. an element of ClK[2].*/

K = bnfinit('x^2+23*47);
ClK = redrepshnf(K);
Cl2 = apply(qfbtohnf,two_torsion(K.disc));
