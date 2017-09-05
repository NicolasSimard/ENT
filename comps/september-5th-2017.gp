K = bnfinit(x^2+47);
pipdata = pipinit(K);
ida = idealprimedec(K,13)[1];
ell = 1;

pol(n) = delkformal('E2,n);

d2l_1E2(pipdata,ell,ida)-dnE2(pipdata,2*ell-1,[ida],"qexp")[1]

qhc = [[1],[2*ell,0]];
print("algo = 0:      ", pnorm(pipdata,qhc));
print("algo = lpsi2:  ", pnorm(pipdata,qhc,"lpsi2"));
print("algo = closure:", pnorm(pipdata,qhc,pol));
print("algo = qexp:   ", pnorm(pipdata,qhc,"qexp"));
print("algo = qexp2:   ", pnorm(pipdata,qhc,"qexp2"));
