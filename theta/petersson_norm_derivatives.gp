/* This script uses he formula for the Petersson norm of theta series in terms
of linear combinations of derivatives of the Eisenstein series E2.*/

E2num(z) = (8*Pi*imag(z))^-1-1/24+suminf(n=1,sigma(n)*exp(2*Pi*I*n*z));

E4num(z) = 1/240+suminf(n=1,sigma(n,3)*exp(2*Pi*I*n*z));

E6num(z) = -1/504+suminf(n=1,sigma(n,5)*exp(2*Pi*I*n*z));

dnE2(n) =
{
    my(v,d);
    v=['E2,'E4,'E6];
    d=[5/6*'E4-2*'E2^2,7/10*'E6-8*'E2*'E4,400*'E4^2-12*'E2*'E6];
    z->subst(subst(subst(diffop('E2,v,d,n),'E2,E2num(z)),'E4,E4num(z)),'E6,E6num(z));
}

dn(P,n) =
{
    my(v,d);
    v=['E2,'E4,'E6];
    d=[5/6*'E4-2*'E2^2,7/10*'E6-8*'E2*'E4,400*'E4^2-12*'E2*'E6];
    diffop(P,v,d,n);
}

dnnum(P,n) =
{
    my(v,d);
    v=['E2,'E4,'E6];
    d=[5/6*'E4-2*'E2^2,7/10*'E6-8*'E2*'E4,400*'E4^2-12*'E2*'E6];
    z->subst(subst(subst(diffop(P,v,d,n),'E2,E2num(z)),'E4,E4num(z)),'E6,E6num(z));
}

pip(pipdata,ell,ida,idb) = 
{
    my(K = pipdata[1]);
    'C_K*idealnorm(K,idb)^(2*ell)*S(pipdata,ell,idealmul(K,ida,idealinv(K,idb)));
}

S(pipdata,ell,ida) =
{
    my(mu, K=pipdata[1], amb = pipdata[3], coords, sqroot, c0);
    coords = bnfisprincipal(K,ida,0);
    for(i=1,#coords,if(coords[i]%2 != 0, return(0)));
    sqroot = vector(#coords,i,coords[i]/2)~;
    c0 = idealinv(K,idealfactorback(K,K.gen,sqroot));
    ac0 = subst(K.zk*bnfisprincipal(K,idealmul(K,idealpow(K,c0,2),ida))[2],'x,K.roots[1]);
    ac0^(2*ell)*sum(i=1,#amb, amb[i][2]^(2*ell)*d2l_1E2(pipdata,ell,idealmul(K,c0,amb[i][1])));
}

\\ Evaluates d^(2*ell-1)E_2 at quadratic ideals
d2l_1E2(pipdata,ell,ida) =
{
    my(mu, K=pipdata[1], i0=1, coords = bnfisprincipal(pipdata[1],ida,0));
    while(pipdata[2][i0][2] != coords, i0 += 1);
    mu = bnfisprincipal(K,idealmul(K,idealinv(K,pipdata[2][i0][1]),ida))[2];
    mu = subst(K.zk*mu,'x,K.roots[1]); \\ ida = mu*pipdata[2][i0][1]
    mu^(-4*ell)*subst(subst(subst(dn('E2,2*ell-1),'E2,pipdata[4][1][i0]),'E4,pipdata[4][2][i0]),'E6,pipdata[4][3][i0]);
}

\\ Petersson inner product init
pipinit(Kraw) =
{
    my(tmp,D=Kraw.disc, K=bnfinit('x^2-D),hK=K.clgp.no);
    my(reps = vector(hK), amb = [], fs, eiseval = vector(3,n,vector(hK)));
    
    fs = reduced_forms(D);
    
    for(i=1,#fs,
        reps[i] = [qfbtohnf(fs[i]),0];
        reps[i][2] = bnfisprincipal(K,reps[i][1],0); \\ Only the coordinates
        tmp = bnfisprincipal(K,idealpow(K,reps[i][1],2));
        if(tmp[1] == 0,
            amb = concat(amb,[[reps[i][1],subst(K.zk*tmp[2],'x,K.roots[1])]]);
        );
    );
    print(reps);
    print(amb);
    
    \\ Evaluate the Eisenstein series at CM points
    print("Evaluating the Eisenstein series...");
    for(i=1,hK,
        tmp = subst(K.zk*reps[i][1],'x,K.roots[1]); \\ tmp = [a,(-b+sqrt(D))/2]
        eiseval[1][i] = tmp[1]^-2*E2num(tmp[2]/tmp[1]);
        eiseval[2][i] = tmp[1]^-4*E4num(tmp[2]/tmp[1]);
        eiseval[3][i] = tmp[1]^-6*E6num(tmp[2]/tmp[1]);
    );
    print("Done.");
    
    return([K,reps,amb,eiseval]);
}