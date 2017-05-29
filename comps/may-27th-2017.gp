/*Script to test the canperiod function.*/

K = bnfinit('x^2+47);
h_K = K.clgp.no;
ida = redrepshnf(K)[1];
[Om, M] = canperiod(K,ida,1);

print("E_4(Om*M)-E_4(ida): ",E(4,Om*M)-E(4,idatolat(K,ida)));
print("\n\nE_6(Om*M)-E_6(ida): ",E(6,Om*M)-E(6,idatolat(K,ida)));
print("\n\ndelta(Om*M)-delta(ida): ",delta(Om*M)-delta(idatolat(K,ida)));

print("\n\nFor E_2: ",factor(algdep((2*Pi*I*Om)^2*E(2,idatolat(K,ida)),2*h_K)));
print("\n\nFor E_6: ",factor(algdep((2*Pi*I*Om)^6*E(6,idatolat(K,ida)),2*h_K)));
print("\n\nFor E_10: ",factor(algdep((2*Pi*I*Om)^10*E(10,idatolat(K,ida)),2*h_K)));
print("\n\nFor delta: ",factor(algdep((2*Pi*I*Om)^12*delta(idatolat(K,ida)),2*h_K)));
