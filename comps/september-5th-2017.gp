K = bnfinit(x^2+47);
pipdata = pipinit(K);
ida = idealprimedec(K,13)[1];
ell = 50;

pol(n) = delkformal('E2,n);

\\d2l_1E2(pipdata,ell,ida)-dnE2(pipdata,2*ell-1,[ida],"qexp")[1];

qhc = [[1],[2*ell,0]];
/*print("algo = 0:      ", N1 = pnorm(pipdata,qhc));
print("algo = lpsi2:  ", N2 = pnorm(pipdata,qhc,"lpsi2"));
print("algo = closure:", N3 = pnorm(pipdata,qhc,pol));*/
print("algo = def:");
N1 = pnorm(pipdata,qhc,"qexp")
/*print("algo = lpsi2:");
N2 = pnorm(pipdata,qhc,"lpsi2")
print("algo = closure:");
N3 = pnorm(pipdata,qhc,pol)*/
print("algo = qexp:");
N4 = pnorm(pipdata,qhc,"qexp")
print("algo = qexpv:");
N5 = pnorm(pipdata,qhc,"qexpv")
