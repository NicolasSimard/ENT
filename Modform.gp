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
    - d(f('q),n) -> d^n(f('q)) (d='q*d/d'q)
    
    *Other functions:
    - jpol(p('q)) -> P(X): P(j('q)) = p('q)
    ");
}

/* To define the Eisenstein series, we follow Zagier's convention:

G_k(s) = 1/2*sum_{(m,n)\neq(0,0)}(mz+n)^-k

E_k(s)=G_k(s)/zeta(k) (leading coefficient is 1)

GG_k(s) = (k-1)!/(2*Pi*I)^kG_k(s) (all non-constant Fourrier coeff are int)
*/

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

\\ Auxiliary funtion to compute delta('q) efficiently
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
    j=j_qexp(k);
    js=vector(k+1,n,j^(k+1-n));
    M=matrix(k+1,k+1,n,m,polcoeff(js[m],n-k-1,'q));
    B=vector(k+1,n,polcoeff(f,n-k-1,'q));
    return(Pol(matsolve(M,B~),'X));
}
{
    addhelp(jpol,"jpol(f): Given a q-expansion f(q) with enough coefficients,
    returns a polynomial g(Y) such that g(j)=f, where j is the j-function.");
}

fd(i:small, pr:small = 15) =
{
    if(i%4 != 0 && i%4 != 3, error("Wrong value for i: has to be cong to 0 or 3 mod 4."));
    my(fis = vector(i+1), tmp, j4); \\ fis[i] = f_{i-1}, so f_i = fis[i+1]
    if(i == 0, return(theta_qexp(pr)), fis[1] = theta_qexp(pr));
    tmp = (theta_qexp(pr)*d(V(4)(E_qexp(10,pr)))/2-10*d(theta_qexp(pr))*V(4)(E_qexp(10,pr)))/V(4)(delta_qexp(pr));
    tmp = (tmp+608*fis[1])/-20;
    if(i == 3, return(tmp), fis[4] = tmp);
    j4 = V(4)(j_qexp(pr));
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

fd2(i:small, pr:small = 20) =
{
    if(i%4 != 0 && i%4 != 3, error("Wrong value for i: has to be cong to 0 or 3 mod 4."));
    my(g1, g1quot, g4, g4quot, j4, j4pow, f0, f3, fi = 0);
    if(i == 0, return(theta_qexp(pr)), f0 = theta_qexp(pr));
    f3 = (theta_qexp(pr)*d(V(4)(E_qexp(10,pr)))/2-10*d(theta_qexp(pr))*V(4)(E_qexp(10,pr)))/V(4)(delta_qexp(pr));
    f3 = (f3+608*f0)/-20;
    if(i == 3, return(f3));
    
    j4 = V(4)(j_qexp(pr));
    g1 = theta1_qexp(pr)*V(4)(E_qexp(4,pr))/V(4)(eta3_qexp(pr))^2/'q;
    g4 = (10*d(g1)*V(4)(E_qexp(10,pr))-3/2*g1*d(V(4)(E_qexp(10,pr))))/V(4)(delta_qexp(pr));
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

gD(D:small, pr:small = 15) =
{
    if(D%4 != 0 && D%4 != 1, error("Wrong value for D: has to be cong to 0 or 1 mod 4."));
    my(gDs = vector(D), tmp, j4); \\ gDs[D] = g_D
    
    gDs[1] = theta1_qexp(pr)*V(4)(E_qexp(4,pr))/V(4)(eta3_qexp(pr))^2/'q;
    if(D == 1, return(gDs[1]));
    
    j4 = V(4)(j_qexp(pr));
    gDs[4] = (10*('q*gDs[1]')*V(4)(E_qexp(10,pr))-3/2*gDs[1]*('q*V(4)(E_qexp(10,pr))'))/V(4)(delta_qexp(pr));
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


