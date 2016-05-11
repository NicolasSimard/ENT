isdisc(D) = return(D%4 <= 1);

control(D) = if(D%4 > 1 || D>=0,error(D," is not a negative discriminant."));

/*Every discriminant uniquely determines an order in a quadratc field. This
function returns its conductor (the index of this order in the maximal order).
*/
conductor(D) =
{
    my(L = factor(D),p);
    control(D);
    p=prod(n=1,length(L~),L[,1][n]^floor(L[,2][n]/2)); \\ D=p^2m, m sq-free
    if(isdisc(D/p^2),return(p),return(p/2));
}

isprimitive(D) = (conductor(D) == 1);

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

sqrt_mod(D,M) =
{
    my(n=1,betas=[]);
    while(n <= floor(M/2),if((n^2-D)%M == 0,betas=concat(betas,[n])); n+=1;);
    return(betas);
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

class_nbr(D) = #primitive_reduced_forms(D);

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

norm_class_nbr(D) = 2*class_nbr(D)/wD(D);

HKclass_nbr(D) =
{
    my(forms = reduced_forms(D));
    sum(i=1,#forms,2/wQ(forms[i]));
}
addhelp(HKclass_nbr,"HKclass_nbr(D): Returns the Hurwitk-Kronecker class number for D, which is defined as the sum over all reduced quandratic forms (not just the primitive ones) of 2/wQ(f) is the order of the stabiliser of the forms f in SL_2(Z).");

genus_nbr(D) =
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

ida(f) = [f[1],(f[2]+sqrt(f[2]^2-4*f[1]*f[3]))/2];
addhelp(ida,"ida(f): Return the ideal attached to the quadratic form f=[a,b,c], i.e."\
"ida(f)=[a,(b+sqrt(D))/2].");

two_torsion(D) =
{
    my(forms, fv);
    forms = primitive_reduced_forms(D);
    fv = [forms[1]];
    for(i=2, length(forms),
        if(forms[i][2] == 0 || forms[i][1] == forms[i][3] || abs(forms[i][2]) == forms[i][1],
        fv = concat(fv,[forms[i]]);
        );
    );
    return(fv);
}

\\ Returns the Chowla-selberg period of discriminant D, as defined in 1-2-3 of
\\ modular forms by Zagier.
CSperiodZag(D) =
{
    control(D);
    return(prod(j=1,abs(D)-1,gamma(j/abs(D))^kronecker(D,j))^(wD(D)/4/class_nbr(D))/sqrt(2*Pi*abs(D)));
}

\\ Returns the Chowla-selberg period of discriminant D, as defined Cohen's book
on Number Theory, volume 2. Gives better results for E2 and E6...
CSperiod(D) =
{
    control(D);
    return(sqrt(prod(j=1,abs(D)-1,gamma(j/abs(D))^kronecker(D,j))^(wD(D)/2/class_nbr(D))/(4*Pi*sqrt(abs(D)))));
}
