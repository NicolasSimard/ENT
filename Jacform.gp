\r modforms.gp

psi101(npr=4) =
{
    my(Dqexp = delta_qexp(npr+1),tau=n->polcoeff(Dqexp,n,'q),Fexp);
    my(c101=vector(4*npr));
    c101[3] = 1;
    for(d=4,4*npr,
        if(d%4 == 0,
            c101[d] = -2*sum(r=1,sqrt(d),if(d-r^2 == 0,0,c101[d-r^2]));
        );
        if(d%4 == 3,
            c101[d] = tau((d+1)/4)-sum(r=2,sqrt(d+1),if(d+1-r^2 == 0,0,r^2*c101[d+1-r^2]));
        )
    );
    Fexp=0;
    for(n=1,npr-1,
        for(r=-floor(2*sqrt(n))+issquare(n),floor(2*sqrt(n))-issquare(n),
            Fexp += c101[4*n-r^2]*'q^n*'z^r)
    );
    Fexp;
}

psi121(npr=4) =
{
    my(Dqexp = delta_qexp(npr+1),tau=n->polcoeff(Dqexp,n,'q),Fexp);
    my(c121=vector(4*npr));
    c121[3] = 1;
    for(d=4,4*npr,
        if(d%4 == 0,
            c121[d] = 12*tau(d/4)-2*sum(r=1,sqrt(d),if(d-r^2 == 0,0,c121[d-r^2]));
        );
        if(d%4 == 3,
            c121[d] = (d+1)/4*tau((d+1)/4)-sum(r=2,sqrt(d+1),if(d+1-r^2 == 0,0,r^2*c121[d+1-r^2]));
        )
    );
    Fexp=0;
    for(n=1,npr-1,
        for(r=-floor(2*sqrt(n))+issquare(n),floor(2*sqrt(n))-issquare(n),
            Fexp += c121[4*n-r^2]*'q^n*'z^r)
    );
    Fexp;
}
