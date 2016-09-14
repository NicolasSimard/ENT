f(K,ida) = {
    my(w,tau);
    w = if(imag(K.roots[1])>0,K.roots[1],conj(K.roots[1]));
    tau = subst(K.zk*idealhnf(K,ida),variable(K),w);
    sqrt(imag(tau[2]/tau[1]))*abs(eta(tau[2]/tau[1],1))^2;
}

eps(K,charcomp) = {
    my(qhcdata=qhcinit(K),reps=redrepshnf(qhcdata[1]),sqchar);
    sqchar=vector(#charcomp,i,2*charcomp[i]);
    prod(n=1,#reps,f(qhcdata[1],reps[n])^qhceval(qhcdata,[sqchar,[0,0]],reps[n]))^-1; 
}
