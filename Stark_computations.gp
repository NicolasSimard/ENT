/*Given an ideal in K, returns the corresponding point in the upper half plane.*/
z_ida(K,ida) = {
    my(tmp=subst(K.zk*idealhnf(K,ida),variable(K),K.roots[1]));
    if(imag(tmp[2]/tmp[1])>0,tmp[2]/tmp[1],tmp[1]/tmp[2]);
}

/*Given a principal ideal, returns a generator.*/
generator(K,ida) = {
    my(tmp=bnfisprincipal(K,ida));
    if(tmp[1] != vector(#tmp[1])~, error("Non-principal ideal."));
    subst(K.zk*tmp[2],variable(K),K.roots[1]);
}

/*F is defined as in Siegel - Lectures in advanced Analytic Number Theory, p.65.*/
F(K,ida) = {
    my(tau=z_ida(K,ida));
    sqrt(imag(tau))*abs(eta(tau,1))^2;
}

/*Defined as in Zagier's chapter in 1-2-3 of Modular Forms.*/
Zag_psi(K,ida) = {
    my(tau=z_ida(K,ida));
    imag(tau)*abs(eta(tau,1))^4;
}

/*rho is defined as in Siegel - Lectures in advanced Analytic Number Theory, p.72.*/
rho(K,ida) = eta(z_ida(K,1),1)^(24*K.clgp.no)/generator(K,idealpow(K,ida,K.clgp.no))^12/eta(z_ida(K,ida),1)^(24*K.clgp.no);

/*logprod is the product of logarithms that appears in the formula for the
Petersson norm of weight one theta series.*/
logprod(K,charcomp) = {
    my(qhcdata=qhcinit(K),reps=redrepshnf(qhcdata[1]),sqchar);
    sqchar=vector(#charcomp,i,2*charcomp[i]);
    prod(n=1,#reps,F(qhcdata[1],reps[n])^qhceval(qhcdata,[sqchar,[0,0]],reps[n]))^-1; 
}
