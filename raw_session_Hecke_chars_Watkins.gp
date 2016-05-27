\p 100
K23=bnfinit(x^2+23);
K23.clgp
K23.zk
p2=idealprimedec(K23,2)[1]
p23=idealpow(K23,p2,3)
bnfisprincipal(K23,p23)
nfbasistoalg(K23,%[2])
mu=(3+sqrt(-23))/2
norm(mu)
z3=exp(2*Pi*I/3);
p2=idealprimedec(K23,2)[1]
C=(2*Pi*I)^2;
idealhnf(K23,p2)
idealhnf(K23,idealpow(K23,p2,2))
tau0=(sqrt(-23)+1)/2; taup2=(sqrt(-23)+3)/4; taup22=(sqrt(-23)+3)/8;
F(i)=C*(E2num(tau0)+z3^i*mu^(2/3)*E2num(taup2)/4+z3^(2*i)*mu^(4/3)*E2num(taup22)/16);