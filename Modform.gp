{
    addhelp(Modform,"This package defines some common modular forms, some operators
    on them and tools to compute basis of certain spaces of modular forms.
    
    In general, a modular form is represented by a q-expansion. The precision
    of the q-expansion is determined by the constant default(seriesprecision)
    which can be changed with the command \ps prec.
    
    *General modular forms:
    - G(k,'q) -> q-expansion; G(k,z) -> G_k(z) (numerical value)
    - E(k,'q) -> q-expansion; E(k,z) -> E_k(z) (numerical value)
    - theta0(x) -> q-expansion
    - theta1(x) -> q-expansion
    - ellj('q) -> q-expansion (PARI built-in); ellj(z) -> j(z) (numerical value)
    - delta('q) -> q-expansion; delta(z) -> delta(z) (numerical value)
    - fd(d) -> q-expansion
    - gD(D) -> q-expansion
    
    *Basis of modular forms:
    - vmbasis(k,N,flag) -> [f1('q),f2('q),...,fd('q)] (q-expansions of Victor Miller basis)
    
    *Operators:
    - U(n,f('q)) -> U_n(f('q)) (U_n operator in level 1)
    - V(n,f('q)) -> V_n(f{'q)) (V_n operator in level 1)
    - pstab(p,f('q)) -> f('q)^[p] (p-stabilisation of f)
    - rcbracket(f('q),k_f,g('q),k_g) -> [f('q),g('q)] (Rankin-Cohen bracket)
    - dop(f('q),n) -> d^n(f('q)) (d='q*d/d'q)
    
    *Other functions:
    - jpol(p('q)) -> P(X): P(j('q)) = p('q)
    ");
}

/*------------------------General modular forms ----------------------*/

G(k,x) =
{
    if(type(x) == "t_POL",
        -bernfrac(k)/2/k
        +subst(Ser(concat([0],vector(default(seriesprecision)-1,n,sigma(n,k-1))),'X),'X,x)
        +O(variable(x)^default(seriesprecision))
    ,
        if(k == 2,
        -1/24+suminf(n=1,sigma(n)*exp(2*Pi*I*n*x)),
        -bernfrac(k)/(2*k)+suminf(n=1,sigma(n,k-1)*exp(2*Pi*I*n*x))
        )
    );
}
{
    addhelp(G,"G(k,x): If x is a polonomial, returns the q-expansion of the
    Eisenstein series G_k(x) (normalized so that the coefficient of 'q is 1).
    If x is anything else, tries to numerically evaluate G_k at that point.");
}

E(k,x) =
{
    if(type(x) == "t_POL",
        1-2*k/bernfrac(k)*subst(Ser(concat([0],vector(default(seriesprecision)-1,n,sigma(n,k-1))),'X),'X,x)
        +O(variable(x)^default(seriesprecision))
    ,
        if(k == 2,
        1-24*suminf(n=1,sigma(n)*exp(2*Pi*I*n*x)),
        1-2*k/bernfrac(k)*suminf(n=1,sigma(n,k-1)*exp(2*Pi*I*n*x))
        )
    );
}
{
    addhelp(E,"E(k,x): If x is a polonomial, returns the q-expansion of the
    Eisenstein series E_k(x) (normalized so that the constant term is 1).
    If x is anything else, tries to numerically evaluate E_k at that point.");
}

E2star(z) = (8*Pi*imag(z))^-1-1/24+suminf(n=1,sigma(n)*exp(2*Pi*I*n*z));
{
    addhelp(E2star,"E2star(z): evaluates the weight 2 Eisenstein series at z.
    E2star is (8*Pi*imag(z))^-1-1/24+O('q).");
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
    ,
        eta(x,1)^24
    );
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
    addhelp(jpol,"jpol(f): Given a q-expansion f(q) with enough coefficients,
    returns a polynomial g(Y) such that g(j)=f, where j is the j-function.");
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
    given, returns the nth iteration of dop.");
}

delkformal(P,n=1) = 
{
    my(v,d);
    v=['E2s,'E4,'E6];
    d=[5/6*'E4-2*'E2s^2,7/10*'E6-8*'E2s*'E4,400/7*'E4^2-12*'E2s*'E6];
    diffop(P,v,d,n);
}
{
    addhelp(delkformal,"delkformal(P,{n=1}): Formal differentiation of an almost holomorphic
    modular form of level 1 represented by a weighted homogeneous polynomial P in
    'E2s, 'E4 and 'E6, where 'E2s is the almost holomorphic weight 2 modular form
    given by E2-3/(Pi*y). The operator delk = q*d/dq - k/(4*Pi*y) preserves this
    ring and the formulas are given in Shimura - Elementary Dirichlet series and
    modular forms. If n is given, returns the nth iteration of del.");
}

pstab(p,f) = f-V(p,U(p,f));
addhelp(pstab,"pstab(p,f): p-stabilisation operator *^[p]=1-U_pV_p on modular forms.");

rcbracket(f,k,g,l) = k*dop(g)*f-l*dop(f)*g;
{
    addhelp(rcbracket,"rcbracket(f,k,g,l): Return the Rankin-Cohen bracket of
    the two modular forms f and g of weigh k and l, respectively.");
}