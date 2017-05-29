/*Script to test the canperiod function.*/

K = bnfinit('x^2+13*7);
h_K = K.clgp.no;
ida = redrepshnf(K)[1];
[Om, M] = canperiod(K,ida,1);

print("E_4(Om*M)-E_4(ida): ",E(4,Om*M)-E(4,idatolat(K,ida)));
print("E_6(Om*M)-E_6(ida): ",E(6,Om*M)-E(6,idatolat(K,ida)));
print("delta(Om*M)-delta(ida): ",delta(Om*M)-delta(idatolat(K,ida)));

print("For E_6: ",factor(algdep((2*Pi*I*Om)^6*E(6,idatolat(K,ida)),2*h_K)));
print("For E_10: ",factor(algdep((2*Pi*I*Om)^10*E(10,idatolat(K,ida)),2*h_K)));
print("For delta: ",factor(algdep((2*Pi*I*Om)^12*delta(idatolat(K,ida)),2*h_K)));
