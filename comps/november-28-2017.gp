deltaprod(D) = {
    my(K = bnfinit('x^2-D), reps = redrepshnf(K));
    prod(i = 1, #reps, idealnorm(K,reps[i])^6*abs(delta(idatolat(K, reps[i]))));
}

CMprod(D, f, k) = {
    my(K = bnfinit('x^2-D), reps = redrepshnf(K));
    prod(i = 1, #reps, idealnorm(K, reps[i])^(k/2)*abs(f(idatolat(K, reps[i]))));
}

gammaprod(D) = {
    prod(j = 1, abs(D), gamma(j/abs(D))^kronecker(D,j))^6/(2 * Pi * abs(D))^(6 * qfbclassno(D));
}

CSperiod2(D) = prod(j = 1, abs(D), gamma(j/abs(D))^kronecker(D,j))^(1 / 2 / qfbclassno(D))/ sqrt(2 * Pi * abs(D));

CSperiodnorm(D) = prod(j = 1, abs(D), gamma(j/abs(D))^kronecker(D,j))^(1 / 2 / qfbclassno(D))/ sqrt(2 * Pi * abs(D)^2);

alpha(D, f, k) = CMprod(D,f,k)/CSperiod2(D)^(k * qfbclassno(D));

alphanorm(D, f, k) = CMprod(D,f,k)/CSperiodnorm(D)^(k*qfbclassno(D));

alpha24(D, f, k) = alpha(D, f, k)^24;

alpha24norm(D, f, k) = alphanorm(D, f, k)^24;

/*
gp > \p 2000
gp > algdep(alpha(-23, x->240*E(4,x), 4),10)
%1 = x - 23375
gp > algdep(alpha(-23, x->504*E(6,x), 6),10)
%2 = x^2 - 14301909271727
gp > algdep(alpha(-23, x->24*E(2,x), 2),10)
%3 = 12167*x^2 - 175561
(13:22) gp > algdep(alpha(-23, x->delta(x), 12),10)
%4 = x - 1

and

gp > \p 2000
gp > algdep(alphanorm(-23, x->24*E(2,x), 2),10)
%1 = x^2 - 2136050687
gp > algdep(alphanorm(-23, x->504*E(6,x), 6),10)
%2 = x^2 - 46397551977112434801878045924234367263
gp > algdep(alphanorm(-23, x->240*E(4,x), 4),10)
%3 = x - 3460338905375
gp > algdep(alphanorm(-23, x->delta(x), 12),10)
%4 = x - 3244150909895248285300369
*/