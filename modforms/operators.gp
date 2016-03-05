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

d =
{
    f -> if(type(f) == "t_SER",
        return('q*f');,
        if(type(f) == "t_CLOSURE",
            return(n -> n*f(n));,
            error("Wrong type for operator d: ",type(f));
        );
    );
}
addhelp(d,"d: Returns the operator d=q*d/dq on modular forms of level 1. Takes modular forms as input\n(represented as a power series or a closure) and returns a modular form represented in the same way.");

add(f,g) = {}

mul(f,g) = {}

bracket(p:small) =
{
    f -> add(f,mul(-1,V(p)(U(p)(f))));
}
