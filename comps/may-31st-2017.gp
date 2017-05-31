/* Script to test if the CM values of the ordinary p-stabilisation of E2 are algebraic. */

K = bnfinit('x^2+23);
h_K = K.clgp.no;
p = 7;
n = 0;
ida = redrepshnf(K)[1];
[Om, M] = canperiod(K,ida,1);

VnE2p(n,p,z) = E(2,p^n*z) - p*E(2,p^(n+1)*z);

\\print("\n\nFor V^",n,"E_2^(",p,"): ",algdep((2*Pi*I)^2*VnE2p(n,p,M),h_K)); \\Not algebraic, it seems
print("For V^",n,"E_2^(",p,"): ",f=algdep((2*Pi*I*Om)^2*VnE2p(n,p,idatolat(K,ida)),h_K));\\ Algebraic
