{
    addhelp(Modform,"This package defines some common modular forms, some operators on them and tools to compute basis of certain spaces of modular forms.
    
    In general, a modular form is represented by a q-expansion. The precision of the q-expansion is determined by the constant default(seriesprecision) which can be changed with the command \ps prec.
    
    *General modular forms:
    - G(k,'q) -> q-expansion; G(k,z) -> G_k(z) in C; G(k,L) -> G_k(L) in C; 
    - E(k,'q) -> q-expansion; E(k,z) -> E_k(z) in C; E(k,L) -> E_k(L) in C;
    -                         G2star(z) -> value in C; G2star(L) -> value in C; 
    - theta0(x) -> q-expansion
    - theta1(x) -> q-expansion
    - ellj('q) -> q-expansion (PARI built-in); ellj(z) -> j(z) (numerical value)
    - delta('q) -> q-expansion; delta(z) -> delta(z) in C; delta(L) -> delta(L) in C;
    - fd(d) -> q-expansion
    - gD(D) -> q-expansion
    
    *Basis of modular forms:
    - vmbasis(k,N,flag) -> [f1('q),f2('q),...,fd('q)] (q-expansions of Victor Miller basis)
    - GktoG4G6(k) -> P(G4,G6) : P(G4,G6)=Gk
    
    *Operators:
    - U(n,f('q)) -> U_n(f('q)) (U_n operator in level 1)
    - V(n,f('q)) -> V_n(f{'q)) (V_n operator in level 1)
    - pstab(p,f('q)) -> f('q)^[p] (p-stabilisation of f)
    - rcbracket(f('q),k_f,g('q),k_g) -> [f('q),g('q)] (Rankin-Cohen bracket)
    - dop(f('q),n) -> d^n(f('q)) (d='q*d/d'q)
    - dopformal(P,n) -> d^n(P) (P in C['E2,'E4,'E6])
    - delkformal(P,n) -> del_k^n(P) (P in C['E2s,'E4,'E6, 'G2s,'G4,'G6])
    
    *Other functions:
    - jpol(p('q)) -> P(X): P(j('q)) = p('q)
    - area([w_1,w_2]) -> area of the lattice Z*w_1+Z*w_2
    ");
}

/*------------------------General modular forms ----------------------*/

/*Note that
elleisnum([w1,w2],k)=(2*Pi*I/w2)^k(1+2/zeta(1-k)*\sum_{n\geq1}sigma_{k-1}(n)q^n)
                    =(2*Pi*I/w2)^k(1+2/zeta(1-k)q+O(q^2))
and
zeta(1-k)=-bernfrac(k)/k.

Note also that G(k,z) is \mathbb{G}_k(z) in Zagier - 1-2-3 of modular forms.

All Eisenstein series were tested in two ways. First, their q-expansion was
compared with know q-expansion up to small precision. Second, the numerical
value of E2star, E4 and E6 at lattices was compared with the CM-values table in Zagier - 1-2-3 of modular forms (after proposition 27). The CM values of delta
were also compared with the entries of the same table. 
*/
G(k,x) =
{
    if(type(x) == "t_POL",
        -bernfrac(k)/2/k + subst(Ser(concat([0],vector(default(seriesprecision)-1,n,sigma(n,k-1))),'X),'X,x) + O(variable(x)^default(seriesprecision))
    , if(type(x) == "t_VEC",
        elleisnum(x,k)/(2*Pi*I)^k/2*zeta(1-k)
    , \\ If x is not a polynomial or a lattice, assume x is complex
        elleisnum([x,1],k)/(2*Pi*I)^k/2*zeta(1-k)
    ));
}
{
    addhelp(G,"G(k,x): If x is a polynomial, returns the q-expansion of the Eisenstein series G_k('q) (normalized so that the coefficient of 'q is 1). This is E_k in Shimura - Elementary Dirichlet series and L-function. If x = [w1,w2] is a Z-basis of a lattice in C, G(k,x) = w2^-k*G_k(w1/w2). Otherwise, assume x is in C and evaluate G_k(x).");
}

E(k,x) = -G(k,x)*2*k/bernfrac(k);
{
    addhelp(E,"E(k,x): If x is a polynomial, returns the q-expansion of the Eisenstein series E_k(x) (normalized so that the constant term is 1). If x = [w1,w2] is a Z-basis of a lattice in C, E(k,x) = w2^-k*E_k(w1/w2). Otherwise, assume x is in C and evaluate E_k(x).");
}

G2star(x) = 
{
    if(type(x) == "t_VEC",
        if(imag(x[1]/x[2])>0,x[2]^-2*G2star(x[1]/x[2]),x[1]^-2*G2star(x[2]/x[1]))
    , \\ If x is not a lattice, assume x is complex
        1/(8*Pi*imag(x))+G(2,x)
    );
}
{
    addhelp(G2star,"G2star(x): If x = [w1,w2] is a Z-basis of a lattice in C, G2star(x) = w2^-2*G_2^*(w1/w2). Otherwise, assume x is in C. G2star(z) is (8*Pi*imag(z))^-1-1/24+O('q).");
}

E2star(x) = -24*G2star(x);
{
    addhelp(E2star,"E2star(x): If x = [w1,w2] is a Z-basis of a lattice in C, E2star(x) = w2^-2*E_2^*(w1/w2). Otherwise, assume x is in C. E2star(z) is -3/(Pi*imag(z))+1+O('q). This is the same E_2^* as in Zagier - 1-2-3 of modular forms.");
}

theta0(x) =
{
    my(pr=default(seriesprecision));
    1+subst(Ser(concat([0],vector(pr-1,n,2*issquare(n))),'X),'X,x)+O(variable(x)^pr);
}
{
    addhelp(theta0,"theta0(x): Returns the q-expansion of the theta series
    (=sum_{n\in\Z} q^(n^2)) of weight 1/2 evaluated at x.");
}

theta1(x) = 
{
    my(d, pr=default(seriesprecision));
    1+subst(Ser(concat([0],vector(pr-1,n,2*issquare(n,&d)*(-1)^d)),'X),'X,x)+O(variable(x)^pr);
}
{
    addhelp(theta1,"theta1(x): Returns the q-expansion of the theta series
    (=sum_{n\in\Z} (-1)^n*q^(n^2)) of weight 1/2 up to precision pr.");
}

\\ Auxiliary function to compute delta('q) efficiently
eta3(x) =
{
    my(d, pr=default(seriesprecision));
    subst(Ser(vector(pr-1,n,if(issquare(1+8*(n-1),&d),(-1)^((-1+d)/2)*d,0)),'X),'X,x)+O(variable(x)^pr);
}

delta(x) =
{
    if(type(x) == "t_POL",
       x*sqr(sqr(sqr(eta3(x))))+O(variable(x)^default(seriesprecision))
    , if(type(x) == "t_VEC",
        if(imag(x[2]/x[1]) > 0, x[1]^-12*eta(x[2]/x[1],1)^24, x[2]^-12*eta(x[1]/x[2],1)^24)
    , \\ If x is not a polynomial or a lattice, assume x is complex
        eta(x,1)^24
    ));
}
addhelp(delta,"delta(x): modular form delta='q*prod(n=1,oo,1-'q^n)^24.");

/*j_qexp(x) = E(4,x)^3/delta_qexp(x);*/

jpol(f) =
{
    my(j,js,M,B,k=-valuation(f,'q));
    if(k <=  0, return(f));
    j=ellj('q);
    js=vector(k+1,n,j^(k+1-n));
    M=matrix(k+1,k+1,n,m,polcoeff(js[m],n-k-1,'q));
    B=vector(k+1,n,polcoeff(f,n-k-1,'q));
    return(Pol(matsolve(M,B~),'X));
}
{
    addhelp(jpol,"jpol(f): Given a q-expansion f(q) with enough coefficients,returns a polynomial g(Y) such that g(j)=f, where j is the j-function.");
}

fd(i:small) =
{
    if(i%4 != 0 && i%4 != 3, error("Wrong value for i: has to be cong to 0 or 3 mod 4."));
    my(fis = vector(i+1), tmp, j4, pr=default(seriesprecision)); \\ fis[i] = f_{i-1}, so f_i = fis[i+1]
    if(i == 0, return(theta0('q)), fis[1] = theta0('q));
    tmp = (theta0('q)*dop(V(4,E(10,'q)))/2-10*dop(theta0('q))*V(4,E(10,'q)))/V(4,delta('q));
    tmp = (tmp+608*fis[1])/-20;
    if(i == 3, return(tmp), fis[4] = tmp);
    j4 = V(4,ellj('q));
    for(j = 4, i,
        if(j%4 == 0 || j%4 == 3,
            tmp = j4*fis[j-4+1];
            for(n = -j + 1, 0,
                if((-n)%4 == 0 || (-n)%4 == 3, 
                    tmp -= polcoeff(tmp,n,'q)*fis[-n+1]
                );
            );
            fis[j+1] = tmp;
        );
    );
    fis[i+1];
}

fd2(i:small) =
{
    if(i%4 != 0 && i%4 != 3, error("Wrong value for i: has to be cong to 0 or 3 mod 4."));
    my(g1, g1quot, g4, g4quot, j4, j4pow, f0, f3, fi = 0, pr=default(seriesprecision));
    if(i == 0, return(theta0('q)), f0 = theta0('q));
    f3 = (theta0('q)*dop(V(4,E(10,'q)))/2-10*dop(theta0('q))*V(4,E(10,'q)))/V(4,delta('q));
    f3 = (f3+608*f0)/-20;
    if(i == 3, return(f3));
    
    j4 = V(4,ellj('q));
    g1 = theta1('q)*V(4,E(4,'q))/V(4,eta3('q))^2/'q;
    g4 = (10*dop(g1)*V(4,E(10,'q))-3/2*g1*dop(V(4,E(10,'q))))/V(4,delta('q));
    g4 = (g4 + (10*j4-21344)*g1)/-20;

    j4pow = 1;
    g1quot = g1/j4;
    g4quot = g4/j4;

    for(n =0, i\4 - (i%4 == 0),
        fi += polcoeff(g4quot,i,'q)*j4pow*f0 + polcoeff(g1quot,i,'q)*j4pow*f3;
        j4pow *= j4;
        g1quot /= j4;
        g4quot /= j4;
    );
    fi += polcoeff(g4quot,i,'q)*j4pow*f0;
}

gD(D:small) =
{
    if(D%4 != 0 && D%4 != 1, error("Wrong value for D: has to be cong to 0 or 1 mod 4."));
    my(gDs = vector(D), tmp, j4, pr=default(seriesprecision)); \\ gDs[D] = g_D
    
    gDs[1] = theta1('q)*V(4,E(4,'q))/V(4,eta3('q))^2/'q;
    if(D == 1, return(gDs[1]));
    
    j4 = V(4,ellj('q));
    gDs[4] = (10*('q*gDs[1]')*V(4,E(10,'q))-3/2*gDs[1]*('q*V(4,E(10,'q))'))/V(4,delta('q));
    gDs[4] = (gDs[4] + (10*j4-21344)*gDs[1])/-20;

    for(j = 5, D,
        if(j%4 == 0 || j%4 == 1,
            tmp = j4*gDs[j-4];
            for(n = -j, -1,
                if((-n)%4 == 0 || (-n)%4 == 1, 
                    tmp -= polcoeff(tmp,n,'q)*gDs[-n]
                );
            );
            gDs[j] = tmp;
        );
    );
    gDs[D];
}

/*------------------------Basis of modular forms ----------------------*/

vmbasis(k,N=0,reduce=1) = 
{
    my(n,e,ls,A,E6,E6_squared,D,Eprod,Dprod,pr=default(seriesprecision));
    if(k%2 == 1,error("The weight ",k," must be even."));
    if(k==0, return([1+O('q^pr)]));

    e=k%12 + 12*(k%12 == 2);
    n=(k-e)/12;

    if(pr <= n,
        ls = vector(n+1);
        ls[1] = 1+ O('q^pr);
        for(i=2,pr,ls[i] = q^i + O('q^pr));
        for(i=pr+1,n+1,ls[i] = O('q^pr));
        return(ls);
    );

    E6 = E(6,'q);
    A = if(e==0,1,E(e,'q)); \\ Slight difference with e=6
    if(polcoeff(A,0,'q) == -1, A=-A);
    D = delta('q);

    if(N>0,E6 = E6*Mod(1,N); A = A*Mod(1,N); D = D*Mod(1,N));

    E6_squared = sqr(E6) + O('q^pr);
    Eprod = E6_squared;
    Dprod = D;

    ls = vector(n+1,i,A);
    for(i=1,n,
        ls[n-i+1] = ls[n-i+1]*Eprod + O('q^pr);
        ls[i+1] = ls[i+1]*Dprod + O('q^pr);
        Eprod = Eprod*E6_squared + O('q^pr);
        Dprod = Dprod*D + O('q^pr);
    );
    \\At this point the basis is upper triangular
    if(reduce == 0,return(ls));
    for(i=1,n,
        for(j=1,i,
            ls[j] = ls[j] - polcoeff(ls[j],i,'q)*ls[i+1];
        );
    );
    return(ls);
}
{
    addhelp(vmbasis,"vmbasis(k,{N=0},{reduce=1}): Returns the Victor Miller
    basis of weight k up to precision pr. If N > 0, returns this basis 
    modulo N. If reduce = 0 (default is 1), the basis is not reduced.");
}

GpoltoEpol(P) = subst(subst(subst(subst(P,'G2s,-'E2s/24),'G2,-'E2/24),'G4,'E4/240),'G6,-'E6/504);

EpoltoGpol(P) = subst(subst(subst(subst(P,'E2s,-24*'G2s),'E2,-24*'E2),'E4,240*'G4),'E6,-504*'G6);

GktoG4G6(k) =
{
    if(k%2 == 1, return(0));
    if(k == 4, return('G4));
    if(k == 6, return('G6));
    6*(k-2)!/(k/2-3)/(k+1)*sum(r=4,k-2,GktoG4G6(r)/(r-2)!*GktoG4G6(k-r)/(k-r-2)!);
}
addhelp(GktoG4G6,"GktoG4G6(k): Express G_k as a polynomial in 'G4 and 'G6.");

EktoE4E6(k) = -2*k/bernfrac(k)*GpoltoEpol(GktoG4G6(k));;
addhelp(EktoE4E6,"EktoE4E6(k): Express E_k as a polynomial in 'E4 and 'E6.");

/*------------------------Operators on modular forms ----------------------*/
U(m,f) =
{
    my(v=variable(f),pr=poldegree(truncate(f),v));
    precision(Ser(vector(floor(pr/m),n,polcoeff(f,(n-1)*m,v)),v),pr);
}
{
    addhelp(U,"U(m,f): Returns the operator U_m on modular forms of level 1.
    Takes a modular form sum_n a_n*q^n as input and returns sum_n a_{mn}*q^n.");
}

V(m,f) =
{
    my(v=variable(f));
    precision(subst(f,v,v^m),default(seriesprecision));
}
{
    addhelp(V,"V(m,f): Returns the operator V_m on modular forms of level 1.
    Takes a modular form sum_n a_n*q^n as input and returns sum_n a_n*q^{mn}.");
}

dop(f,t=1) =
{
    my(v=variable(f),pr=poldegree(truncate(f),v));
    Ser(concat([0],vector(pr-1,n,n^t*polcoeff(f,n,v))),v)+O(v^pr);
}
{
    addhelp(dop,"dop(f,{t=1}): Operator d=q*d/dq on modular forms of level 1.
    Takes a modular form sum_n a_n*q^n and returns sum_n n^t*a_n*q^n.");
}

dopformal(P,n=1) =
{
    my(v,d);
    v=['E2,'E4,'E6];
    d=[('E2^2-'E4)/12,('E2*'E4-'E6)/3,('E2*'E6-'E4^2)/2];
    diffop(P,v,d,n);
}
{
    addhelp(dopformal,"dopformal(P,{n=1}): Formal differentiation of a quasimodular
    form of level 1 represented by a weighted homogeneous polynomial P in
    'E2, 'E4 and 'E6. The operator dop = q*d/dq preserves this ring and the
    formulas are given in Zagier - 1,2,3 of modular forms (prop 15). If n is
    given, returns the nth iteration of dopformal. Note that E_k = 1 + O(q).");
}

delkformal(P,n=1) = 
{
    my(v,d);
    v=['G2s,'G4,'G6,'E2s,'E4,'E6];
    dG = [5/6*'G4-2*'G2s^2,7/10*'G6-8*'G2s*'G4,400/7*'G4^2-12*'G2s*'G6];
    dE = [('E2s^2-'E4)/12,('E4*'E2s-'E6)/3,('E6*'E2s-'E4^2)/2];
    diffop(P,v,concat(dG,dE),n);
}
{
    addhelp(delkformal,"delkformal(P,{n=1}): Formal differentiation of an almost holomorphic
    modular form of level 1 represented by a weighted homogeneous polynomial P in
    'E2s, 'E4 and 'E6, where 'E2s is the almost holomorphic weight 2 modular form
    given by E2-3/(Pi*y), or a weighted homogeneous polynomial P in 'G2s, 'G4
    and 'G6. The operator delk = q*d/dq - k/(4*Pi*y) preserves this ring and
    the formulas are given in Shimura - Elementary Dirichlet series and
    modular forms. If n is given, returns the nth iteration of delkformal.
    Note that G_k = -B_k/2k + O(q) and E_k = 1 + O(q).");
}

pstab(p,f) = f-V(p,U(p,f));
addhelp(pstab,"pstab(p,f): p-stabilisation operator *^[p]=1-U_pV_p on modular forms.");

rcbracket(f,k,g,l) = k*dop(g)*f-l*dop(f)*g;
{
    addhelp(rcbracket,"rcbracket(f,k,g,l): Return the Rankin-Cohen bracket of
    the two modular forms f and g of weigh k and l, respectively.");
}

/*--------------- Other functions --------------------------------------*/
area(L) = matdet([real(L[1]),imag(L[1]); real(L[2]),imag(L[2])]);
{
    addhelp(area,"area(L): return the area of the Z-lattice L. If L = [w1,w2] is a 2-component vector of linearly independent complex numbers, returns the area of the lattice spanned by those numbers in C.");
}