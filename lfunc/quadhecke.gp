qhlinit(K) =
{
    my(hK = K.clgp.no, reps = redrepshnf(K));    
    return([qhcinit(K),reps,matrix(hK,3,i,j,E(2*j,idatolat(K,reps[i])))]);
}
addhelp(qhlinit,"qhlinit(K): initialize data to evaluate the L-function of a Hecke character of an imaginary quadratic field. This is a vector of the form [qhcinit(K),reps,eiseval].");

/* Tested: qhlfun(qhlinit(bnfinit('x^2+23)),[[i],[2,0]],2) = values in Watkin's paper

and

prod(i=0,2,qhlfun(qhlinit(bnfinit('x^2+23)),[[i],[2,0]],2))=
Pi^2/5*qhlfun(qhlinit(bnfinit('x^2+23)),[[0],[6,0]],4)

as in Watkin's paper.*/
qhlfun(qhldata,qhc,m) =
{
    \\ for the moment, we only accept T of the form [2*ell,0]
    if(qhc[2][2] != 0 || qhc[2][1]%2 == 1, error("Wrong values for infinity type: ",qhc[2],"."));
    
    my(t = qhc[2][1], K=qhldata[1][1], hK=K.clgp.no, reps = qhldata[2]);
    
    \\ Make sure qhc is not the trivial class character
    if(t == 0 && vector(#qhc[1],i,qhc[1][i]%K.clgp[2][i]) == vector(#qhc[1]),
            error("L-function of K for trivial character has a pole at s=1.")
    );
    
    \\ If qhc is a trivial character, make sure m=1
    if(t == 0 && m != 1, error("Not implemented yet at s != 1 for class characters."));
    
    \\ Make sure m is in the correct range for non class characters
    if(t != 0 && (2*m-t < 2 || t-m < 0), error("Not implemented for those values of m."));
    
    \\ Compute the value!
    if(t == 0,
        -Pi/6/sqrt(abs(K.disc))*sum(i=1,#reps,qhceval(qhldata[1],qhc,reps[i])*log(idealnorm(K,reps[i])^6*abs(delta(idatolat(K,reps[i])))));
        ,
        my(S, mf);
    
        mf = if(2*m-t == 2, delkformal('E2,t-m), delkformal(EktoE4E6(2*m-t),t-m));
        S = sum(i=1,hK,qhceval(qhldata[1],qhc,qhldata[2][i])*substvec(mf,['E2,'E4,'E6],qhldata[3][i,]));
        I^t*(2*Pi)^m/gamma(m)*sqrt(abs(K.disc))^(t-m)*S;
    );
}
addhelp(qhlfun,"qhlfun(qhldata,qhc,m): Evaluates the Hecke L-function attached to the Hecke character [c,T] (where c are the components and T is the infinity type) at the point m (for the moment, we need t/2+1 <= m <= t, or t = 0 and m = 1). qhldata is the data returned by qhlinit.");
