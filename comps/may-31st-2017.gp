/* Script to test if the CM values of the ordinary p-stabilisation of E2 are algebraic. */

K = bnfinit('x^2+23);
h_K = K.clgp.no;
p = 13;
idp = idealprimedec(K,p)[1];
ida = redrepshnf(K)[2];
[Om, M] = canperiod(K,ida,1);

idpmnida(n) = idealmul(K,idealpow(K,idp,-n),ida);
VnE2p(n) = E(2,idatouhp(K,idpmnida(n))) - p*E(2,idatouhp(K,idpmnida(n+1)));

\\print("\n\nFor V^",n,"E_2^(",p,"): ",algdep((2*Pi*I)^2*VnE2p(n,p,M),h_K)); \\Not algebraic, it seems
\\print("For V^",n,"E_2^(",p,"): ",f=algdep((2*Pi*I*Om)^2*VnE2p(n,p,idatolat(K,ida)),h_K)); \\ Algebraic
print("p=",p,". Is inert: ",(#idealprimedec(K,p) == 1));
for(m=0,10,f=algdep((2*Pi*I*canperiod(K,idpmnida(n)))^2*VnE2p(m),2*h_K); print("degree:",poldegree(f),". f=",factor(f),". \n Valuation:"))