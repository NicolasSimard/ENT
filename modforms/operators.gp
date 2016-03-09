U(p:small) = 
{
    f -> if(type(f) == "t_SER",
        my(pr=poldegree(truncate(f),'q)+1);
        return(Ser(vector(floor((pr-1)/p),n,polcoeff(f,(n-1)*p,'q)),'q));,
        if(type(f) == "t_CLOSURE",
            return(n -> f(p*n));,
            error("Wrong type for operator U_",p,":",type(f));
        );
    );
}
addhelp(U,"U(p): Returns the operator U_p on modular forms of level 1. Takes modular forms as input\n(represented as a power series or a closure) and returns a modular form represented in the same way.");

V(p:small) =
{
    f -> if(type(f) == "t_SER",
        return(substpol(f,'q,'q^p));,
        if(type(f) == "t_CLOSURE",
            return(n -> if(n%p != 0,0,f(n/p)));,
            error("Wrong type for operator V_",p,":",type(f));
        );
    );
}
addhelp(V,"V(p): Returns the operator V_p on modular forms of level 1. Takes modular forms as input\n(represented as a power series or a closure) and returns a modular form represented in the same way.");

d(f,t=0) =
{
    if(type(f) == "t_SER",
        if(t < 0 || type(t) != "t_INT", error("Not implemented yet!"));
        if(t == 0, return('q*f'), return(d(f,t-1))),
        if(type(f) == "t_CLOSURE",
            if(t < 0 && f(0) != 0, error("Has to be a cusp form"));
            return(n -> if(n == 0, 0, f(n)*n^t));,
            error("Wrong type for operator d: ",type(f));
        );
    );
}
addhelp(d,"d(f,{t=0}): Operator d=q*d/dq on modular forms of level 1. Takes modular forms as input\n(represented as a power series or a closure) and returns a modular form represented in the same way. If the optionnal parameter t is >0, applies the d operator t times. If t < 0, applies the inverse operator (formal integration), which makes sense for p-adic cusp forms.");

mfadd(f,g) = 
{
    if(type(f) == "t_SER" && type(g) == "t_SER",
        return(f+g);,
        if(type(f) == "t_CLOSURE" && type(g) == "t_CLOSURE",
            return(n -> f(n) + g(n));,
            error("Wrong type:",type(f)," and ",type(g),".");
        );
    );
}
addhelp(mfadd,"mfadd(f,g): Returns the sum of the modular forms f and g. Returns the same type as the input (either a power series or a closure).");

mfmul(f,g) = 
{
    if(type(f) == "t_SER" && type(g) == "t_SER", return(f*g));
    if(type(f) == "t_CLOSURE" && type(g) == "t_CLOSURE",return(n -> sum(k=0,n,f(k)*g(n-k))));
    if(type(f) == "t_SER" || type(g) == "t_SER", return(f*g));
    if(type(f) == "t_CLOSURE", return(n -> f(n)*g));
    if(type(g) == "t_CLOSURE", return(n -> f*g(n)));
    error("Wrong type:",type(f)," and ",type(g),".");
}
addhelp(mfmul,"mfmul(f,g): Returns the product of two modular forms or the product of a scalar and a modular form. Returns the same type as the input (either a power series or a closure).");

bracket(p:small) = f -> mfadd(f,mfmul(-1,V(p)(U(p)(f))));
addhelp(bracket,"bracket(p): Returns the bracket operator *^[p]=1-U_pV_p on modular forms.");

rankincohen(f,k:small,g,l:small) = mfadd(mfmul(k,mfmul(f,d(g))),mfmul(-l,mfmul(d(f),g)));
addhelp(rankincohen,"rankincohen(f,k,g,l): Return the Rankin-Cohen brancket of the two modular forms f and g of weigh k and l, respectively.");
