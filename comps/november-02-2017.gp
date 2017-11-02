f1_A(K, ida, ell, idb) = {
    my(qexp = Vec(bintheta(K,ida,ell)), L = idatolat(K,idb));
    theval(K, ida, ell, idatolat(K,idb))/E(2, idatolat(K,idealinv(K,ida)))^ell;
}

f2_A(K, ida, ell, idb) = {
    my(qexp = Vec(bintheta(K,ida,ell)), L = idatolat(K,idb));
    theval(K, ida, ell, idatolat(K,idb))/theval(K, 1, ell, idatolat(K,ida));
}