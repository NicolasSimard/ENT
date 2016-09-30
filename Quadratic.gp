{
    addhelp(Quadratic,"This package defines a few functions that deal with
    quadratic forms and quadratic fields. A quadratic form is represented by
    a 3-component vector [a,b,c]. In PARI, they are represented by objects of
    type t_QFI created as Qfb(a,b,c). 
    
    *general functions
    - control(D) -> 1 or error
    - conductor(D) -> f: D=f^2D_0
    - liftSL2Z(M,N) -> M_0: M_0 = M mod N
    - sqrt_mod(D,M) -> [beta]: beta[i]^2 = D mod M 
    
    *Positive definite binary quadratic forms
    - right_act(f,M) -> f([x,y]*M)
    - reduced_forms(D) -> [f]
    - primitive_reduced_forms(D) -> [f]
    - tau(f) -> z
    - reduced_roots(D) -> [tau(f)]
    - primitive_reduced_roots -> [tau(f)]
    - wD(D) -> w_D
    - wQ(f) -> w(Q)
    - qfbclassno(D) -> h_D
    - qfbhclassno(d) -> H(d) (d>0)
    - genusno(D) -> #genus
    - two_torsion(D) -> [amb] = ClK[2]
    - discofclassno(n) -> [D_1,...,D_k]: h(D_i) = n (n <= 10)
    
    *Heegner forms
    - Heegner_form(f,N,beta) -> g: g~f and N|a(g)
    - Heegner_forms(D,N,{beta=[]}) -> Q_{D,N,beta}
    - Heegner_points(D,N,{beta=[]}) -> roots of Q_{D,N,beta}
    
    *Conversion
    - qfbtohnf(f) -> [a,*;0,*]
    - idatoqfb(K,ida) -> [a,b,c]
    
    *Arithmetic
    - CSperiod(D) -> Omega_D
    - CSperiodCoh(D) -> Omega_D
    
    *Imaginary quadratic fields
    - redrepshnf(K) -> [a_1,...,a_h]
    - parirepshnf(K) -> [a_1,...,a_h]
    
    *Hecke characters of Imaginary quadratic fields
    - qhcinit(K) -> [K,[r_i]]
    - qhchars(K,T=[0,0]) -> [psi_1,...,psi_h_K]
    - qhceval(qhcdata,psi,ida) -> psi(ida)
    ");
}

/*------------------------------General functions----------------------------*/

isdisc(D) = return(D%4 <= 1);

control(D) = if(D%4 > 1 || D>=0,error(D," is not a negative discriminant."));

/*Every discriminant uniquely determines an order in a quadratc field. This
function returns its conductor (the index of this order in the maximal order).*/
conductor(D) =
{
    my(L = factor(D),p);
    control(D);
    p=prod(n=1,length(L~),L[,1][n]^floor(L[,2][n]/2)); \\ D=p^2m, m sq-free
    if(isdisc(D/p^2),return(p),return(p/2));
}

right_act(f,M) =
{
    return([f[1]*M[1,1]^2+f[2]*M[1,1]*M[2,1]+f[3]*M[2,1]^2,
            2*f[1]*M[1,1]*M[1,2]+f[2]*(M[1,1]*M[2,2]+M[1,2]*M[2,1])+2*f[3]*M[2,1]*M[2,2],
            f[1]*M[1,2]^2+f[2]*M[1,2]*M[2,2]+f[3]*M[2,2]^2]
    );
}

/* Lift a matrix M in SL2(Z/NZ) to a matrix in SL2(Z).*/
liftSL2Z(M,N) =
{
    my(MZ = lift(M));
    if(gcd(gcd(MZ[2,1],MZ[2,2]),N) != 1,error(M," does not belong to SL2(Z/",N,"Z)."));
    if(N<=1,return(MZ));
    my(mysol,mya,myb,myc=MZ[2,1],myd=MZ[2,2],myu=1,myv=1,myt,myg,myp);
    myg=gcd(myc,myd);
    if(myg != 1,
        myp = factor(myg)[,1];
        if(myc != 0,
            for(n=1,length(myp~),if(myc%myp[n] == 0 && myp[n] != -1,myu*=myp[n],myv*=myp[n]));
            myt=lift(chinese(Mod(1,myu),Mod(0,myv)));
            myd+=myt*N;,
            for(n=1,length(myp~),if(myd%myp[n] == 0 && myp[n] != -1,myu*=myp[n],myv*=myp[n]));
            myt=lift(chinese(Mod(1,myu),Mod(0,myv)));
            myc+=myt*N;
        );
    );
    \\ We now have a matrix mod N with coprime c and d (in Z).
    myc=lift(myc); myd=lift(myd);
    mysol=matsolvemod([myd,-myc;1,0;0,1],[0,N,N]~,[1,MZ[1,1],MZ[1,2]]~)~;
    mya=mysol[1]; myb=mysol[2];
    return([mya,myb;myc,myd]);
}

sqrt_mod(D,M) =
{
    my(n=1,betas=[]);
    while(n <= floor(M/2),if((n^2-D)%M == 0,betas=concat(betas,[n])); n+=1;);
    return(betas);
}

/*-----------------------------Heegner points--------------------------------*/

/* Given N>=0, a discriminant D such that D is a square mod 4N, an integer mod 2N
beta such that beta^2=D mod 4N and a factorisation m1*m2 of gcd(N,beta,(beta^2-D)/4N),
one has an explicit bijection between the sets Q^0_{D,N,beta}/Gamma_0(N) and
Q^0_D/Gamma(1), where Q^0_D is the set of primitive binary quadratic forms of
discriminant D and

Q^0_{D,N,,beta}={[a,b,c] in Q^0_D|a=0 mod N, b=beta mod 2N}.

The bijection is explicit and depends of m1 and m2 (see Gross-Kohnen-Zagier, Heegner
points and derivatives of L-series II). This function takes as input a quadratic form,
an integer N and a beta as above and returns a form in Q^0_{D,N,beta} maping to it under
the bijection. For the moment, the function assumes m1=m2=1, otherwise returns an error.
*/
Heegner_form(f,beta,N,m1=1,m2=1) =
{
    my(sol,m,a,b,c,D,M,r,s,t,u,T);
    a=f[1]; b=f[2]; c=f[3];
    D=b^2-4*a*c;
    m=gcd(gcd(N,beta),(beta^2-D)/4/N);
    if(m != 1, error("Function not implemented for m>1."));
    M=[a,(b-beta)/2;(b+beta)/2,c];
    sol=matsolvemod(M,N,0,1)[2];
    if(sol[,1] != [0,0]~, r=sol[1,1]; t=sol[2,1];, r=sol[1,2]; t=sol[2,2];);
    sol=matsolvemod(Mat(lift([r,-t])),N,1)~;
    u=sol[1]; s=sol[2];
    T=liftSL2Z([r,s;t,u],N);
    return(right_act(f,T));
}

Heegner_forms(D,N,beta=[]) =
{
    control(D);
    if(type(beta) != "t_INT", beta = sqrt_mod(D,4*N), beta=[beta]);
    if(length(beta) == 0, error("The pair (",D,",",N,") does not satisfy the Heegner hypothesis."));
    if((beta[1]^2-D)%(4*N) != 0, error("The integer ",beta[1]," is not a root of ",D," mod 4*",N));
    my(forms = primitive_reduced_forms(D));
    return(vector(length(forms),n,Heegner_form(forms[n],beta[1],N)));
}

Heegner_points(D,N,beta=[]) =
{
    my(Hforms = Heegner_forms(D,N,beta));
    return(vector(length(Hforms),n,tau(Hforms[n])));
}

/*---------------------Binary quadratic forms--------------------------------*/

reduced_forms(D) =
{
    my(b0 = D%2, fv,a,c,zv,n);

    control(D);

    fv = [[1,b0,(b0^2-D)/4]];
    forstep(b = b0, floor(sqrt(-D/3)), 2,
        zv=divisors((b^2-D)/4);
        n=length(zv);
        forstep(j=(n+1)\2, 2, -1,
            a=zv[j]; if (a < b, break);
            c=zv[n-j+1];
            fv = concat(fv,[[a,b,c]]);
            if(b && a != b && a != c, fv = concat(fv,[[a,-b,c]]));
        );
    );
    return(fv);
}

primitive_reduced_forms(D) =
{
    my(forms = reduced_forms(D), prim = []);
    for(i=1, length(forms),
        if(gcd(gcd(forms[i][1],forms[i][2]),forms[i][3]) == 1,
           prim = concat(prim,[forms[i]])
        );
    );
    return(prim);
}

tau(f) = (-f[2]+sqrt(f[2]^2-4*f[1]*f[3]))/2/f[1];

reduced_roots(D) =
{
    my(forms = reduced_forms(D),taus = []);
    for(i=1, length(forms),taus = concat(taus,[tau(forms[i])]));
    return(taus);
}

primitive_reduced_roots(D) =
{
    my(forms=primitive_reduced_forms(D),taus=[]);
    for(i=1, #forms,taus = concat(taus,[tau(forms[i])]));
    return(taus);
}

wD(D) =
{
    if(D%4 > 2, error("Not a discriminant."));
    if(D == -3, return(6),
    if(D == -4, return(4),
                return(2)));
}

wQ(f) =
{
    if(f[1] == f[2] && f[2] == f[3], return(6),
    if(f[1] == f[3] && f[2] == 0,    return(4),
                                     return(2)));
}

norm_class_nbr(D) = 2*qfbclassno(D)/wD(D);

genusno(D) =
{
    my(r,n,mu);
    if(D>=0 || D%4 > 1, return(-1));
    r = omega(-D) - (D%2 == 0);
    if(D%4 == 1, return(2^(r-1)));
    n = floor(-D/4);
    if(n%4 == 3,mu=r);
    if(n%4 == 1 || n%4 == 2 || n%8 == 4,mu=r+1);
    if(n%8 == 0,mu=r+2);
    return(2^(mu-1));
}

two_torsion(D) =
{
    my(forms, fv=[]);
    forms = primitive_reduced_forms(D);
    for(i=1, length(forms),
        if(forms[i][2] == 0 || forms[i][1] == forms[i][3] || abs(forms[i][2]) == forms[i][1],
        fv = concat(fv,[forms[i]]);
        );
    );
    return(fv);
}
addhelp(twotorsion,"two_torsion(D): Return representatives of the two-torsion of the class group of discriminant D. The number"\
"number of such classes is equal to the number of genera for the discriminant D.");

discofclassno(n) =
{
    if(n > 11, error("Not implemented yet. Too lazy!"));
    my(discs = []);
    for(D=3,15667,
        if(isfundamental(-D) && qfbclassno(-D) == n,
            discs = concat(discs,[-D])
        )
    );
    return(discs);
}
{
    addhelp(discofclassno,"discofclassno(n): Return all fundamental discriminants
    with class number n. The problem is solved only for n <= 100 by Watkins,
    and implemented for n <= 11 for the moment.");
}

/*-----------------------------Conversion------------------------------------*/

qfbtohnf(f) =
{
    if(type(f) == "t_QFI", f=Vec(f));
    my(D=f[2]^2-4*f[1]*f[3], K=nfinit('w^2-D), t = nfalgtobasis(K,(-f[2]+'w)/2));
    [f[1],t[1];0,t[2]];
}
{
    addhelp(qfbtohnf,"qfbtohnf(f): Return the Hermite normal form of the ideal corresponding to f in the integral basis of nfinit(x^2-D).");
}

idatoqfb(K,ida) = {
    my(Zbasis = subst(K.zk*idealhnf(K,ida),variable(K),K.roots[1]));
    Vec(round(('X*Zbasis[1]-Zbasis[2])*conj('X*Zbasis[1]-Zbasis[2])/idealnorm(K,ida)));
}
{
    addhelp(idatoqfb,"idatoqfb(K,ida): Given an ideal ida=[a,b] in a quadratic field K, return the corresponding binary quadratic form N(ax-by)/N(ida).");
}

/*-------------------------Imaginary quadratic fields------------------------*/

/* Return the Minkowski bound of a number field K (not necessarily quadratic).*/
Minkowski(K) = my(n=K.r1+2*K.r2); n!/(n^n)*(4/Pi)^K.r2*sqrt(abs(K.disc));

redrepshnf(K) =
{
    my(ClK=K.clgp, reps = []);
    forvec(e=vector(#ClK.cyc,i,[0,ClK.cyc[i]-1]),
        reps=concat(reps,[idealred(K,idealfactorback(K,ClK.gen,e))]);
    );
    reps;
}

parirepshnf(K) =
{
    my(ClK=K.clgp, reps = []);
    forvec(e=vector(#ClK.cyc,i,[0,ClK.cyc[i]-1]),
        reps=concat(reps,[idealfactorback(K,ClK.gen,e)]);
    );
    reps;
}

/*-------------Hecke characters of Imaginary quadratic fields----------------*/
qhcinit(K) =
{
    my(mu, ClK = K.clgp, psi0=vector(#ClK.cyc));
    for(i=1,#ClK.cyc,
        mu = bnfisprincipal(K,idealpow(K,ClK.gen[i],ClK.cyc[i]))[2];
        psi0[i] = subst(K.zk*mu,variable(K),K.roots[1])^(1/ClK.cyc[i]);
    );
    [K,psi0];
}
{
    addhelp(qhcinit,"qhcinit(K): Given an imaginary quadratic field K, computes the principalisations y_i of the generators of ClK, i.e. if g_i is a generator of ClK of order o_i, then g_i^o_i=y_i O_K. Returns [K,[r_i]], where r_i=y_i^(1/o_i).");
}

qhchars(K,T=[0,0]) =
{
    my(qccs=[]);
    forvec(e=vector(#K.clgp.cyc,i,[0,K.clgp.cyc[i]-1]),
        qccs = concat(qccs,[[e,T]]);
    );
    qccs;
}
{
    addhelp(qhchars,"qhchars(K): Given an imaginary quadratic field K, returns all Hecke characters of K of infinity type T, represented as vectors [c,T], where c is the component vector in terms of an elementary divisor decomposition of ClK.");
}

qhceval(qhcdata,qhc,ida) =
{
    my(K=qhcdata[1],ClK = K.clgp, decomp = bnfisprincipal(K,ida), mu);
    if(#qhc[1] != #ClK.cyc, error("Invalid component vector."));
    mu = subst(K.zk*decomp[2],variable(K),K.roots[1]);
    mu^qhc[2][1]*conj(mu^qhc[2][2])*prod(i=1,#ClK.cyc,(qhcdata[2][i]^qhc[2][1]*conj(qhcdata[2][i])^qhc[2][2]*exp(2*Pi*I*qhc[1][i]/ClK.cyc[i]))^decomp[1][i]);
}
{
    addhelp(qhceval,"qhceval(qhdata,qhc,ida): Any Hecke character [c,T] of K evaluated at a generator g_i of the class group has value psi(g_i) = r_i^T[1]*conj(r_i)^T[2]*zeta_i^c_i, where r_i is as above and zeta_i = exp(2*Pi*I/o_i) and 0 <=c_i<o_i. Then if ida = mu*g_1^e_1*...*g_d^e_d, one can determine psi(ida) from the above information. qhdata is returned by qhcinit(K).");
}

/*-----------------------------Arithmetic------------------------------------*/
/*Returns the Chowla-selberg period of discriminant D, as defined in 1-2-3 of
modular forms by Zagier, but with a 4 instead of 2 in the factor.*/
CSperiod(D) =
{
    control(D);
    return(prod(j=1,abs(D)-1,gamma(j/abs(D))^kronecker(D,j))^(wD(D)/4/qfbclassno(D))/sqrt(4*Pi*abs(D)));
}

/* Returns the Chowla-selberg period of discriminant D, as defined Cohen's book
on Number Theory, volume 2. Gives better results for E2 and E6...*/
CSperiodCoh(D) =
{
    control(D);
    return(sqrt(prod(j=1,abs(D)-1,gamma(j/abs(D))^kronecker(D,j))^(wD(D)/2/qfbclassno(D))/(4*Pi*sqrt(abs(D)))));
}