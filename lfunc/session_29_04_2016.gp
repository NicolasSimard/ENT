K23 = bnfinit(x^2+23);
intBasis = K23.zk;
p2 = idealprimedec(K23,2)[1]; \\ A generator of the class group
lambCoord = bnfisprincipal(K23,idealpow(K23,p2,3))[2]; \\ principalisation
lamb = subst(intBasis*lambCoord,x,sqrt(-23));

\\ We take as class reps Z_K, p2, p2^2.
\\ We first find a Z-basis of each of them.

ida1=subst(intBasis*idealpow(K23,p2,0,1),x,sqrt(-23));\\Z_K
ida2=subst(intBasis*idealpow(K23,p2,1,1),x,sqrt(-23));\\p2
ida3=subst(intBasis*idealpow(K23,p2,2,1),x,sqrt(-23));\\p2^2

char(n,t) = [[ida1,1],[ida2,exp(2*Pi*I*n/3)*lamb^(t/3)],[ida3,(exp(2*Pi*I*n/3)*lamb^(t/3))^2]];
