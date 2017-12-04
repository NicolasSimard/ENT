read("peterssonnorm.gp"); \\ For lsym2level1data

default(realprecision, 250); \\ Set real precision
default(seriesprecision, 10000); \\ number of coeffs in the q-exps of delta and E(k)

delta = eta(q)^24;
E(k) = 1 - 2 * k / bernfrac(k) * Ser(concat([0], vector(default(seriesprecision) - 1, n, sigma(n, k - 1))), 'q) + O('q^default(seriesprecision));

k = 12;
coeffs = Vec(delta); \\ q + O(q^2)*/

/*k = 16;
coeffs = Vec(delta * E(4)); \\ q + O(q^2)*/

/*k = 18;
coeffs = Vec(delta * E(6)); \\ q + O(q^2)*/

/*k = 20;
coeffs = Vec(delta * E(8)); \\ q + O(q^2)*/

/*k = 22;
coeffs = Vec(delta * E(10)); \\ q + O(q^2)*/

/*k = 26; \\ Not in their paper, but dim S_26(SL_2(Z)) = 1
coeffs = Vec(delta * E(14)); \\ q + O(q^2)*/

print(">>> Initializing the L-function. Can take a few seconds/minutes...");
/* The argument [(3 * k - 2) / 2, (k - 2) / 2, 0] in lfuninit defines a range of
real values containing k and 2*k-2, the two values we are interested in...*/
L = lfuninit(lfuncreate(lsym2level1data(coeffs, k)), [(3 * k - 2) / 2, (k - 2) / 2, 0]);
c(k) = Pi^(4 - 2 * k) * 2^(2 * k - 2 * k + 2) * (k / 2 - 1)! * (2 * k - 3)! / (k - 1)!;

print(">>> Running algdep. Polynomial saved as f.");
print(f = algdep(c(k) * lfun(L, 2 * k - 2) / lfun(L, k), 10, round(0.8 * default(realprecision))));

{if(poldegree(f) == 1,
    print(">>> f has degree 1, so the quantity is rational and we factor it:");
    print(factor(-polcoeff(f, 0) / polcoeff(f, 1))));
}