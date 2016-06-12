\\ For idoneal numbers
{
    my(idoneal=[]);
    for(n=1,1848,
        if(isfundamental(-4*n) && qfbclassno(-4*n) == genusno(-4*n),idoneal=concat(idoneal,[-4*n]))
    );
    for(i=1,#idoneal,
        print(idoneal[i],":");
        dataD = pipinit(bnfinit(x^2-idoneal[i]));
        periodD = CSperiod(idoneal[i]);
        for(ell=1,5,print(Vec(algdep(pip(dataD,ell,1,1)/periodD^(4*ell),3))))
    )
}
/*
\\For class number one
{
my(dataD,periodD);
for(D=5,163,
    if(isfundamental(-D) && qfbclassno(-D)==1,
        dataD = pipinit(bnfinit(x^2+D));
        periodD = CSperiod(-D);
        for(ell=1,5,print(Vec(algdep(pip(dataD,ell,1,1)/periodD^(4*ell),3))))
    )
)
}

\\For class number 2
{
my(dataD,periodD);
for(D=5,427,
    if(isfundamental(-D) && qfbclassno(-D) == 2,
        dataD = pipinit(bnfinit(x^2+D));
        periodD = CSperiod(-D);
        print(D,":");
        for(ell=1,5,print(Vec(algdep(pip(dataD,ell,1,1)/periodD^(4*ell),3))))
    )
)
}*/