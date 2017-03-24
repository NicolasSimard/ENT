/*L-function of Hecke characters.*/

/* Initial parameters. */
D = -47;
K = bnfinit(x^2-D);
G = K.clgp;

ell = 0;
T = [2*ell,0];
comps = vector(#G.cyc,i,random(G.cyc[i]));
hchar = [comps,T];

n = 11 ;

binthetaqhc(K,hchar,'q)