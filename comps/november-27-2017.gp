p(ell) = my(char = qhchars(K,[2*ell,0])); prod(i = 1, #reps, pnorm(pipdata, char[i]));

p2(ell) = my(char = qhchars(K,[2*ell,0])); prod(i = 1, #reps, lfun(lfuncreate(qhcLdata(K, char[i])),ell+1));

\\ algdep(p(ell)/pnorm(pipdata,[[0],[2*ell*#reps,0]]),10) \\Rational
algdep(p2(ell)/lfun(lfuncreate(qhcLdata(K,  [[0],[2*ell*#reps,0]])),ell*#reps + 1)/Pi^(#reps-1),10) \\Rational
