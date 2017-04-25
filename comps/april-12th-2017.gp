D = -47;
K = bnfinit(x^2-D);
ida = idealprimedec(K,7)[1];
phi_ida = deltaquot(K,ida);