phi(K,ida) = {
    my(L = idatolat(K, ida));
    sqrt(idealnorm(K, ida)) * abs(L[1]^-1 * eta(L[2] / L[1], 1)^2);
}

kappa_psi(K, c) = {
    my(reps = redrepshnf(K));
    prod(i = 1, #reps, phi(K, reps[i])^(-qhceval(qhcinit(K), [c, [0, 0]], reps[i])^2));
}

test(D, deg) = {
    my(c, K = bnfinit('x^2-D));
    print("Class group: ", K.clgp);
    if(!deg, deg = 2 * K.clgp.no);
    
    while(1,
        print(">>> Input a class char: ", end = "");
        c = input();
        if(c,
            if(qhcistrivial(K, 2 * [c,[0, 0]]),
                print("!!! Genus character")
            ,
                print(pol = algdep(kappa_psi(K, c), deg)));
        ,
            break()));
}

quadhilbertoverq(D) = polredbest(polcompositum('x^2-D, quadhilbert(D))[1]);