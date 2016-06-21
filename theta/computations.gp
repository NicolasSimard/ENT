factortex(N) = 
{
    my(M,str="", tmp);
    M=factor(N);
    for(i=1,#M[,1],
        if(M[i,1] == -1,
            tmp = "-",
            if(M[i,2] == 1,
                tmp = if(i == #M[,1],Str(M[i,1]),concat(Str(M[i,1]),"\\cdot")),
                tmp = concat(concat(concat(Str(M[i,1]),"^{"),Str(M[i,2])),"}");
            )
        );
        str=concat(str,tmp);
    );
    return(concat(concat("$",str),"$"));
}

\\ For idoneal numbers (Note: -35 is not idoneal, but clgp has exponent 2...)
idoneal=[];
for(n=1,1848,if(isfundamental(-4*n) && qfbclassno(-4*n) == genusno(-4*n) && qfbclassno(-4*n) > 2, idoneal=concat(idoneal,[-4*n])));

table(discs, maxell, out) =
{
print("Number of discs: ",#discs);
my(comps = [],dataD,periodD,tmp,D);
for(i=1,#discs,
    D = discs[i];
    dataD = pipinit(bnfinit(x^2-D));
    periodD = CSperiod(D);
    comps = concat(comps,[concat(concat("$",Str(D)),"$")]);
    print(D);
    for(ell=1,maxell,
        tmp = Vec(algdep(pip(dataD,ell,1,1)/periodD^(4*ell),2));
        comps = concat(comps,[if(#tmp == 2, factortex(tmp[2]/tmp[1]), "X")]);
    )
);
/*
writetex(out,
    matrix(#comps/(maxell+1),maxell+1,D,ell,comps[(D-1)*(maxell+1)+ell])
);*/
for(i=1,#comps,write(out,comps[i]));
}
