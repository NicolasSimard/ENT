f_A(K, ida, ell, idb) = {
    my(qexp = Vec(bintheta(K,ida,ell)), L = idatolat(K,idb));
    L[1]^-(2*ell+1)*sum(n = 1, #qexp, qexp[n] * exp(2 * Pi * I * n * L[2]/L[1]))*E(2, idatolat(K,ida))^ell;
}