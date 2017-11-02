{addhelp(Thetafunc,"This package defines the binary theta functions attached to imaginary quadratic fields. More precisely, it defines functions to compute q-expansions.

*compute q-expansions
- bintheta(K,X,ell,var='q,prec=0,flag=0) -> power series");
}

bintheta(K, X, ell, prec) = {
    if(type(X) == "t_VEC" && #X == 2, return(binthetaqhc(K, X, prec)));
    return(binthetaida(K, X, ell, prec));
}
{addhelp(bintheta,"bintheta(K, X, ell, {prec = auto}): Return the q-expansion of the binary theta series attached to K, X and ell. Here, X can be a Hecke character of K (i.e. X = [[c_1,...,c_d],[2*ell,0]]), in which case ell is ignored and the function returns the q-expansion of the associated newform, or X can be an ideal of K, in which case the corresponding q-expansion is returned.")
}

binthetaeval(K, X, ell, L) = {
    my(qexp = bintheta(K, X, ell), v, coeffs);
    v = valuation(qexp, variable(qexp));
    coeffs = Vec(qexp);
    L[1]^-(2 * ell+1) * sum(n = v, #coeffs, coeffs[n - v + 1] * exp(2 * Pi * I * n * L[2]/L[1]));
}
{addhelp(binthetaeval,"binthetaeval(K, X, ell, L): Evaluate the theta function bintheta(K, X, ell) at the lattice L.")}

/*binthetaqhc1(K,qhc,var='q,flag) = {
    my(L=ideallist(K,default(seriesprecision)), data = qhcinit(K));
    Ser(concat([qhcistrivial(K,qhc)*K.clgp.no/K.tu[1]],vector(#L,n,sum(N=1,#L[n],qhceval(data,qhc,L[n][N])))),var);
}*/

binthetaida(K, ida, ell, prec) = {
    my(abc, L, M, reps, qexp);
    prec = if(prec == 0, default(seriesprecision));
    L = idatolat(K,ida);
    abc = Vec(round(('X * L[1] + L[2]) * conj('X * L[1] + L[2]) / idealnorm(K, ida)));
    M = [abc[1],abc[2]/2;abc[2]/2,abc[3]];
    reps = qfminim(M, prec - 1, , 2)[3];
    /* The following algorithm works, but is much slower...*/
    /*qexp = (ell == 0) + O('q^prec); \\ The constant term
    for(i = 1, #reps[1, ],
        qexp += (reps[1,i]*L[1] + reps[2,i]*L[2])^(2*ell) * 'q^qfeval(M,reps[,i]));
    return(qexp);*/
    coeffs = vector(prec);
    coeffs[1] = ell == 0; \\constant term is 1 iff ell = 0
    for(i = 1, #reps[1, ],
        coeffs[qfeval(M,reps[,i]) + 1] += 2*(reps[1,i]*L[1] + reps[2,i]*L[2])^(2*ell));
    return(Ser(coeffs,'q));
};
{addhelp(binthetaida,"binthetaida(K, ida, ell, prec): Return the q-expansion of the binary theta series attached to K, ida and ell.");
}

binthetaqhc(K, qhc, prec) = {
    my(ClK = redrepshnf(K), data = qhcinit(K), tmp);
    tmp = sum(i = 1, #ClK,
        qhceval(data, qhc, ClK[i])^-1 * binthetaida(K, ClK[i], qhc[2][1]/2, prec))/K.tu[1];
    if(qhc[2] == [0,0] && !qhcistrivial(K,qhc),
        variable(tmp) * Ser(Vec(tmp)[2..-1], variable(tmp), prec)
    ,
        tmp
    )
}
{addhelp(binthetaqhc,"binthetaqhc(K, qhc, prec): Return the q-expansion of the binary theta series attached to K and qhc. Here, qhc is a Hecke character of K (i.e. X = [[c_1,...,c_d],[2*ell,0]]).");
}