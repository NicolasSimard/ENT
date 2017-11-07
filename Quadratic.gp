{
    addhelp(Quadratic,"This package defines a few functions that deal with
    quadratic forms and quadratic fields. A quadratic form is represented by
    a 3-component vector [a,b,c]. In PARI, they are represented by objects of
    type t_QFI created as Qfb(a,b,c). 
    
    *general functions
    - control(D) -> 1 or error
    - conductor(D) -> f: D=f^2D_0
    - factordisc(D) -> [d_1,...,d_k]
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
    - idatouhp(K,ida) -> tau_ida
    - idatolat(K,ida) -> [w1,w2]
    
    *Arithmetic
    - CSperiod(D) -> Omega_D
    - CSperiodCoh(D) -> Omega_D
    - canperiod(K,ida,flag) -> Omega
    
    *Imaginary quadratic fields
    - redrepshnf(K) -> [a_1,...,a_h]
    - parirepshnf(K) -> [a_1,...,a_h]
    - complexgen(K,ida) -> alpha in C s.t. (alpha) = ida
    - idaindex(K,reps,ida) -> [mu,i] s.t. ida = mu*reps[i]
    
    *Hecke characters of imaginary quadratic fields
    - qhcinit(K) -> [K,[r_i]]
    - qhchars(K,T=[0,0]) -> [psi_1,...,psi_h_K]
    - qhceval(qhcdata,psi,ida) -> psi(ida)
    - qhcistrivial(K,qhc) -> 0 or 1
	
	*L-functions attached to imaginary quadratic fields
	- qhcLdata(K,qhc,{eps}) -> Ldata = [a,as,Vga,k,N,eps]
    - lsym2data(K,qhc) -> Lsym2data = [a,as,Vga,k,N,eps]
    ");
}

/*------------------------------General functions----------------------------*/
quadsetup(D, flag) = {
    K = bnfinit('x^2-D);
    reps = redrepshnf(K);
    if(flag >= 1, pipdata = pipinit(K));
    if(flag >= 2, Om_K = CSperiod(K.disc));
    print("Class group: ",K.clgp);
}
addhelp(quadsetup,"quadsetup(D, {flag}): Define K = Q(sqrt(D)) and reps = redrepshnf(K). If flag == 1, define pipdata = pipinit(K) and if flag == 2, define also Om_K = CSperiod(K.disc). This function also print the structure of the class group. ");

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

factordisc(D) =
{
    if(!isfundamental(D) || D>0,error("D has to be a field discriminant."));
    my(v,L);
    v = factor(-D);
    L = vector(#v[,1],i,if(v[i,1] != 2,(-1)^((v[i,1]-1)/2)*v[i,1], 0));
    if(D%2 == 0, L[1] = D/prod(i=2,#L,L[i]));
    L;
}
addhelp(factordisc,"factordisc(D): Return the factorisation of D into prime discriminants.");

right_act(f,M) =
{
    return([f[1]*M[1,1]^2+f[2]*M[1,1]*M[2,1]+f[3]*M[2,1]^2,
            2*f[1]*M[1,1]*M[1,2]+f[2]*(M[1,1]*M[2,2]+M[1,2]*M[2,1])+2*f[3]*M[2,1]*M[2,2],
            f[1]*M[1,2]^2+f[2]*M[1,2]*M[2,2]+f[3]*M[2,2]^2]
    );
}
{
    addhelp(right_act,"right_act(f,M): right action of the matrix M in M_2(R)
    on the binary quadratic form f.");
}

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
addhelp(liftSL2Z,"liftSL2Z(M,N): lift the matrix M in SL_2(Z/NZ) to SL_2(Z).");

sqrt_mod(D,M) =
{
    my(n=1,betas=[]);
    while(n <= floor(M/2),if((n^2-D)%M == 0,betas=concat(betas,[n])); n+=1;);
    return(betas);
}
addhelp(sqrt_mod,"sqrt_mod(D,M): return all x mod M such that x^2=D mod M.");

/*-----------------------------Heegner points--------------------------------*/
/* Given N>=0, a discriminant D such that D is a square mod 4N, an integer mod 2N
beta such that beta^2=D mod 4N and a factorisation m1*m2 of gcd(N,beta,(beta^2-D)/4N),
one has an explicit bijection between the sets Q^0_{D,N,beta}/Gamma_0(N) and
Q^0_D/Gamma(1), where Q^0_D is the set of primitive binary quadratic forms of
discriminant D and

Q^0_{D,N,beta}={[a,b,c] in Q^0_D|a=0 mod N, b=beta mod 2N}.

The bijection is explicit and depends on m1 and m2 (see Gross-Kohnen-Zagier, Heegner
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
{
    addhelp(Heegner_form,"Heegner_form(f,beta,N,{m1=1,m2=1}): Given N>=0, a
    discriminant D such that D is a square mod 4N, an integer beta mod 2N such
    that beta^2=D mod 4N and a factorisation m1*m2 of gcd(N,beta,(beta^2-D)/4N),
    one has an explicit bijection between the sets Q^0_{D,N,beta}/Gamma_0(N) and
    Q^0_D/Gamma(1), where Q^0_D is the set of primitive binary quadratic forms of
    discriminant D and Q^0_{D,N,beta}={[a,b,c] in Q^0_D|a=0 mod N, b=beta mod 2N}.
    The bijection depends on m1 and m2. Given f and beta, N, m1 and m2 as above,
    returns a form in Q^0_{D,N,beta} mapping to f under the bijection. For the
    moment, the function assumes m1=m2=1 (otherwise returns an error).");
}

Heegner_forms(D,N,beta=[]) =
{
    control(D);
    if(type(beta) != "t_INT", beta = sqrt_mod(D,4*N), beta=[beta]);
    if(length(beta) == 0, error("The pair (",D,",",N,") does not satisfy the Heegner hypothesis."));
    if((beta[1]^2-D)%(4*N) != 0, error("The integer ",beta[1]," is not a root of ",D," mod 4*",N));
    my(forms = primitive_reduced_forms(D));
    return(vector(#forms,n,Heegner_form(forms[n],beta[1],N)));
}
{
    addhelp(Heegner_forms,"Heegner_forms(D,N,{beta=[]}): Return a set of
    representatives of Q^0_{D,N,beta}/Gamma_0(N) (i.e. a set of Heegner forms).
    If beta is not given, we choose one at random. See Heegner_form help for
    more details on Heegner forms.");
}

Heegner_points(D,N,beta=[]) =
{
    my(Hforms = Heegner_forms(D,N,beta));
    return(vector(length(Hforms),n,tau(Hforms[n])));
}
{
    addhelp(Heegner_points,"Heegner_points(D,N,{beta=[]}): Return a set of
    complex points in the complex upper half plane corresponding to a set of
    representatives of Q^0_{D,N,beta}/Gamma_0(N) (i.e. a set of Heegner forms).
    If beta is not given, we choose one at random. See Heegner_form help for
    more details on Heegner forms.");
}

/*---------------------Binary quadratic forms--------------------------------*/
vectoQfb(f) = Qfb(f[1],f[2],f[3]);
addhelp(vectoQfb,"vectoQfb(f): return Qfb(f[1],f[2],f[3]).");

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
{
    addhelp(genusno,"genusno(D): return the number of genus in the class group
    of discriminant D.");
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
{
    addhelp(twotorsion,"two_torsion(D): Return representatives of the two-torsion
    of the class group of discriminant D. The number of such classes is equal to
    the number of genera for the discriminant D.");
}

discofclassno(n) =
{
    if(n > 100, error("The discriminants of class number >100 are not known yet."));
    local(liste);
    read("discs_of_class_no1-100.gp");
    liste[n];
}
{
    addhelp(discofclassno,"discofclassno(n): Return all fundamental discriminants
    with class number n. The problem is solved only for n <= 100 by Watkins.");
}

/*-----------------------------Conversion------------------------------------*/
qfbtohnf(f) =
{
    if(type(f) == "t_QFI", f=Vec(f));
    my(D=f[2]^2-4*f[1]*f[3], K=nfinit('w^2-D), t = nfalgtobasis(K,(-f[2]+'w)/2));
    [f[1],t[1];0,t[2]];
}
{
    addhelp(qfbtohnf,"qfbtohnf(f): Return the Hermite normal form of the ideal
    corresponding to f in the integral basis of nfinit(x^2-D).");
}

idatoqfb(K,ida) =
{
    my(Zbasis = subst(K.zk*idealhnf(K,ida),variable(K),K.roots[1]));
    Vec(round(('X*Zbasis[1]-Zbasis[2])*conj('X*Zbasis[1]-Zbasis[2])/idealnorm(K,ida)));
}
{
    addhelp(idatoqfb,"idatoqfb(K,ida): Given an ideal ida=[a,b] in a quadratic
    field K, return the corresponding binary quadratic form N(ax-by)/N(ida).");
}

idatouhp(K,ida) = {
    my(tmp=subst(K.zk*idealhnf(K,ida),variable(K),K.roots[1]));
    if(imag(tmp[2]/tmp[1])>0,tmp[2]/tmp[1],tmp[1]/tmp[2]);
}
addhelp(idatouhp,"idatouhp(K,ida): Given an ideal ida=[a,b] in a quadratic field K, return the corresponding point b/a or a/b in the upper-half plane.");

idatolat(K,ida) = {
    my(w = if(imag(K.roots[1])>0,K.roots[1],conj(K.roots[1])), L);
    L = subst(K.zk*idealhnf(K,ida),variable(K),w); \\[a,(-b+sqrt(D))/2]
    return(if(imag(L[2]/L[1])<0,L*[0,1;1,0],L));
}
addhelp(idatolat,"idatolat(K,ida): return [w1,w2] such that ida = Z*w1+Z*w2, with w2/w1 in the upper half plane.");

/*-------------------------Imaginary quadratic fields------------------------*/
minkbound(nf) = my(n=nf.r1+2*nf.r2); n!/(n^n)*(4/Pi)^nf.r2*sqrt(abs(nf.disc));
addhelp(minkbound,"minkbound(nf): return the Minkowski bound of a number field nf.");

redrepshnf(K) = apply(ida->idealred(K,ida),parirepshnf(K));
addhelp(redrepshnf,"redrepshnf(K): return a list of reduced ideals generating the ideal class group.");

parirepshnf(K) = {
    my(ClK=K.clgp, reps = []);
    forvec(e=vector(#ClK.cyc,i,[0,ClK.cyc[i]-1]),
        reps=concat(reps,[idealfactorback(K,ClK.gen,e)]));
    reps
}
addhelp(parirepshnf,"parirepshnf(K): return a list of ideals generating the ideal class group.");

complexgen(K,ida) = {
    my(tmp=bnfisprincipal(K,ida));
    if(tmp[1] != vector(#tmp[1])~, error("Non-principal ideal."));
    subst(K.zk*tmp[2],variable(K),K.roots[1])
}
addhelp(complexgen,"complexgen(K,ida): return mu in C such that ida = (mu).");

idarep(K,reps,ida) = {
    my(i,mu);
    i = select(rep->bnfisprincipal(K,idealmul(K,ida,idealinv(K,rep)),0) == vector(#K.clgp.cyc)~, reps, 1)[1];
    mu = complexgen(K,idealmul(K,ida,idealinv(K,reps[i])));
    [mu,i]
}
addhelp(idarep,"idarep(K,reps,ida): return [mu,i] such that ida = mu*reps[i], where mu is a complex number.");

quaddata(D) =
{
    my(data=[[D,factordisc(D)]],K,v,reps);
    K=bnfinit(x^2-D);
    reps = primitive_reduced_forms(D);
    
    \\ Class group part: [h_K, [elementary div], [gens of cyclic comp]]
    data = concat(data,[K.clgp]);
    for(i=1,#data[2][3],data[2][3][i] = Vec(qfbred(vectoQfb(idatoqfb(K,data[2][3][i])))));
    
    \\ Genus theory part: [genus no, [2-torsion Cl_K], [Cl_K^2]]
    v = [genusno(D), two_torsion(D)];
    v = concat(v,[Set(vector(#reps,i,Vec(qfbred(qfbpowraw(vectoQfb(reps[i]),2)))))]);
    data = concat(data,[v]);
    
    \\ Class field theory: [Hilbert class polynomial]
    data = concat(data,[[quadhilbert(D)]]);
}
{
    addhelp(quaddata,"quaddata(D): return the data attached to the imaginary
    quadratic field of discriminant D. The data is of the form
    [[disc],[clgp],[genus theory],[cft]], where
    - [disc] = [D,prime disc facto]
    - [clgp] = [h_D,[cyc],[gen]]
    - [genus theory] = [genus no, Cl_K[2], Cl_k^2]
    - [cft] = [Hlibert class polynomial]");
}

quaddatacheck(data) =
{
    my(v=[]);
    v = concat(v,[data[1][1] == prod(i=1,#data[1][2],data[1][2][i])]);
    v = concat(v,[data[3][1] == 2^(#data[1][2]-1)]);
    v = concat(v,[#data[3][2] == data[2][1]/#data[3][3]]);
    if(prod(i=1,#v,v[i]),1,v);
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
    addhelp(qhcinit,"qhcinit(K): Given an imaginary quadratic field K, computes
    the principalisations y_i of the generators of ClK, i.e. if g_i is a generator
    of ClK of order o_i, then g_i^o_i=y_i O_K. Returns [K,[r_i]],
    where r_i=y_i^(1/o_i).");
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
    addhelp(qhchars,"qhchars(K): Given an imaginary quadratic field K, returns
    all Hecke characters of K of infinity type T, represented as vectors [c,T],
    where c is the component vector in terms of an elementary divisor
    decomposition of ClK.");
}

qhceval(qhcdata,qhc,ida) =
{
    my(K=qhcdata[1],ClK = K.clgp, decomp = bnfisprincipal(K,ida), mu);
    if(#qhc[1] != #ClK.cyc, error("Invalid component vector."));
    mu = subst(K.zk*decomp[2],variable(K),K.roots[1]);
    mu^qhc[2][1]*conj(mu^qhc[2][2])*prod(i=1,#ClK.cyc,(qhcdata[2][i]^qhc[2][1]*conj(qhcdata[2][i])^qhc[2][2]*exp(2*Pi*I*qhc[1][i]/ClK.cyc[i]))^decomp[1][i]);
}
{
    addhelp(qhceval,"qhceval(qhdata,qhc,ida): Any Hecke character [c,T] of K
    evaluated at a generator g_i of the class group has value
    psi(g_i) = r_i^T[1]*conj(r_i)^T[2]*zeta_i^c_i, where r_i is as above and
    zeta_i = exp(2*Pi*I/o_i) and 0 <=c_i<o_i. Then if
    ida = mu*g_1^e_1*...*g_d^e_d, one can determine psi(ida) from the above
    information. qhdata is returned by qhcinit(K).");
}

qhcLdata(K,qhc,eps) =
{
	my(a,Ldata);
	
	a = (n -> Vec(bintheta(K,qhc,1,'q,n)));
	Ldata = [a,1,[0,1],qhc[2][1] + qhc[2][2] + 1,abs(K.disc),'X];
	
	\\ Compute the root number is eps = 0. Otherwise set it to eps
	if(eps == 0, Ldata[6] = lfunrootres(lfuncreate(Ldata))[3], Ldata[6] = eps);
	
	return(Ldata);
}
{
	addhelp(qhcLdata,"qhcLdata(K,qhc,{eps=0}): Return the data [a,astar,Vga,k,N,eps] for the L-function of the Hecke character qhc of K. If eps is not given (i.e. is 0), it will be automatically computed. This can be couputationnally costly and a good guess for eps is  1.");
}

lsym2data(K,qhc) =
{
    if(qhc[2][2] != 0, error("Infinity type must be of the form [2*ell,0]"));
	
    my(ell, k, N, chi_D, anf, a);
    ell = qhc[2][1]/2;
    
    k = 2*ell+1;
    N = abs(K.disc);
    chi_D = (n -> kronecker(K.disc,n));
    
    a = (n-> my(an = Vec(bintheta(K,qhc,1,'q,n)),ap2); 
        direuler(p=2,n,
            if(N%p == 0,
                1/(1-p^(k-1)*X)
            ,
                ap2 = an[p]^2-chi_D(p)*p^(k-1);
                1/(1 - chi_D(p)*ap2*X + chi_D(p)*ap2*p^(k-1)*X^2 - p^(3*(k-1))*X^3))));

    return([a,1,[0,1,2-k],2*k-1,N^2,1]);

}
{
	addhelp(lsym2data,"lsym2data(K,qhc,{eps=0}): Return the data [a,astar,Vga,k,N,eps] for the symmetric square L-function of the theta series attached to the Hecke character qhc = [comp,T] of the number field K.");
}

qhcistrivial(K,qhc) = {
    if(qhc[2] != [0,0], return(0));
    vector(#qhc[1],i,qhc[1][i] % K.clgp.cyc[i]) == vector(#qhc[1]);
}
{
    addhelp(qhcistrivial,"qhcistrivial(K,qhc): Return 1 if the Hecke character is trivial and 0 otherwise.");
}

/*-----------------------------Arithmetic------------------------------------*/
/*Returns the Chowla-selberg period of discriminant D, as defined in 1-2-3 of
modular forms by Zagier, but with a 4 instead of 2 in the factor.*/
CSperiod(D) =
{
    control(D);
    return(prod(j=1,abs(D)-1,gamma(j/abs(D))^kronecker(D,j))^(wD(D)/4/qfbclassno(D))/sqrt(4*Pi*abs(D)));
}
{
    addhelp(CSperiod,"CSperiod(D): return the Chowla-selberg period of
    discriminant D, as defined in 1-2-3 of modular forms by Zagier, but with a
    4 instead of 2 in the factor:
    prod(j=1,|D|-1,gamma(j/|D|)^\chi_D(j))^(wD(D)/4/h_D)/sqrt(4*Pi*|D|)")
}

/* Returns the Chowla-selberg period of discriminant D, as defined Cohen's book
on Number Theory, volume 2. Gives better results for E2 and E6...*/
CSperiodCoh(D) =
{
    control(D);
    return(sqrt(prod(j=1,abs(D)-1,gamma(j/abs(D))^kronecker(D,j))^(wD(D)/2/qfbclassno(D))/(4*Pi*sqrt(abs(D)))));
}
{
    addhelp(CSperiodCoh,"CSperiodCoh(D): return the Chowla-selberg period of
    discriminant D, as defined Cohen's book on Number Theory, volume 2:
    sqrt(prod(j=1,|D|-1,gamma(j/|D|)^\chi_D(j))^(wD(D)/2/h_D)/(4*Pi*sqrt(|D|)))")
}

canperiod(K,ida,flag) =
{
    my(E, idb, j, M, f, t, dist, rootsinK, tau, alpha, Om);
    j = ellj(idatouhp(K,ida));
    E = if(j == 1728,
            ellinit([1,0]),
        if(j == 0,
            ellinit([0,1]),
        \\ else (i.e. j != 0,1728)
            ellinit([1,0,0,-36/(j-1728),-1/(j-1728)])
        ));
    M = ellperiods(E); \\ C/M is defined over the Hilbert class field of K
    
    t = varhigher("t",variable(K));
    f = subst(algdep(M[1]/M[2],2),'x,t);
    rootsinK = nfroots(K,f);
    dist = apply((r)->subst(K.zk*nfalgtobasis(K,r),variable(K),K.roots[1])-M[1]/M[2],rootsinK);
    tau = if(abs(dist[1]) < abs(dist[2]),rootsinK[1], rootsinK[2]); \\ M[1]/M[2] = tau
    idb = idealhnf(K,1,tau);
    alpha = bnfisprincipal(K,idealmul(K,ida,idealinv(K,idb)))[2];
    Om = subst(K.zk*alpha,variable(K), K.roots[1])/M[2];
    
    if(flag,[Om,M],Om);
}
addhelp(canperiod,"canperiod(K,ida,{flag = 0}): Return a somewhat canonical period Om such that ida = Om*M for some lattice M in C such that the elliptic curve C/M is defined over the Hilbert class field of K. This period is well-defined up to a unit in O_K. If flag != 0, return M also.");
