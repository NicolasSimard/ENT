a(n) = vector(n,i,1);
astar = 1;
Vga = [0];
w = 1;
cond = 1;
eps = 1;
res = lfunrootres(lfuncreate([a,astar,Vga,w,cond,eps,0]))[1];
Ldata = [a,astar,Vga,w,cond,eps,res];

printf("Estimated error: %.2f%%", (1.0-abs(lfuncheckfeq(lfuncreate(Ldata)))/default(realbitprecision))*100);
