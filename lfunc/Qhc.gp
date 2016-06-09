qhcinit(K) =
{
    my(mu, ClK = K.clgp, psi0=vector(#ClK.cyc));
    for(i=1,#ClK.cyc,
        mu = bnfisprincipal(K,idealpow(K,ClK.gen[i],ClK.cyc[i]))[2];
        psi0[i] = subst(K.zk*mu,variable(K),K.roots[1]);
    );
    [K,psi0];
}
{
    addhelp(qhcinit,"qhcinit(K): Given an imaginary quadratic field K, computes 
    the principalisations y_i of the generators of ClK, i.e. if g_i is a
    generator of ClK of order o_i, then g_i^o_i=y_i O_K. Returns
    [K,[r_i]], where r_i=y_i^(1/o_i).");
}

qcchars(K) =
{
    my(qccs=[]);
    forvec(e=vector(#K.clgp.cyc,i,[0,K.clgp.cyc[i]-1]),
        qccs = concat(qccs,[[e,[0,0]]]);
    );
    qccs;
}
{
    addhelp(qcchars,"qcchars(K): Given an imaginary quadratic field K, returns
    all characters of the class group of K, represented as vectors [c,T],
    where c is the component vector in terms of an elementary divisor
    decomposition of ClK and T=[0,0] is the infinity type of the character.");
}

qhceval(qhcdata,qhc,ida) =
{
    my(K=qhcdata[1],ClK = K.clgp, decomp = bnfisprincipal(K,ida));
    if(#qhc[1] != #ClK.cyc, error("Invalid component vector."));
    mu = subst(K.zk*decomp[2],variable(K),K.roots[1]);
    mu^qhc[2][1]*conj(mu^qhc[2][2])*prod(i=1,#ClK.cyc,(qhcdata[2][i]^qhc[2][1]*conj(qhcdata[2][i])^qhc[2][2]*exp(2*Pi*I*qhc[1][i]/ClK.cyc[i]))^decomp[1][i]);
}
{
    addhelp(qhchar,"qhchar(qhdata,qhc,ida): Any Hecke character [c,T] of K
    evaluated at a generator g_i of the class group has value
    psi(g_i) = r_i^T[1]*conj(r_i)^T[2]*zeta_i^c_i, where r_i is as above and
    zeta_i = exp(2*Pi*I/o_i) and 0 <=c_i<o_i. Then if
    ida = mu*g_1^e_1*...*g_d^e_d, one can determine psi(ida) from the above
    information. qhdata is returned by qhcinit(K).");
}