/*This script is meant to compute the q-expansions of the theta series
attached to a quadratic field, up to a precision. It implements two functions:

1) reduced_forms(D): return a list of triples representing the reduced forms
of discriminant D.

2) theta_qexps(D, k, N): return a list of the coefficients of the theta series
attached to the quadratic field of discriminant D and with parameter k. Note
that the weight of these theta series is not k.
*/

reduced_forms(D) = {
    local(b0 = D%2, fv,a,c,zv);

    if(D >= 0 || D%4 > 1, return([]));

    fv = [[1,b0,(b0^2-D)/4]];
    forstep(b = b0, floor(sqrt(-D/3)), 2,
        zv=divisors((b^2-D)/4);
        n=length(zv);
        forstep(j=(n+1)\2, 2, -1,
            a=zv[j]; if (a < b, break);
            c=zv[n-j+1];
            if(gcd(gcd(a,b),c) != 1, next);
            fv = concat(fv,[[a,b,c]]);
            if(b && a != b && a != c, fv = concat(fv,[[a,-b,c]]));
        );
    );
    fv
}

theta_qexps(D, k, N)  = {
    local(Bound, forms, qexps, f, n, qexp, alpha, beta);

    if(k == 0, return(0));

    forms = reduced_forms(D);
    qexps = [];

    for(i =1, length(forms),
        f = forms[i];

        \\ Initialise the q-expansion with zeroes
        qexp = [];
        for(j = 1, N-1, qexp = concat(qexp,[0]));

        \\ The generators of the ideal coresponding to f
        alpha = f[1];
        beta = (-f[2]+sqrt(D))/2;

        Bound = floor(sqrt(2*(N-1))/sqrt(f[1]+f[3]-sqrt((f[1]-f[3])^2+f[2]^2)));
        for(x = -Bound, Bound,
            for(y = -Bound, Bound,
                n = f[1]*x^2 + f[2]*x*y + f[3]*y^2;
                if(n < N && n != 0,
                    qexp[n] += (alpha*x - beta*y)^(2*k);
                );
            );
        );

        \\ Note that qexp[n] = c_n, so there is no constant term
        qexps = concat(qexps,[qexp]);
    );
    qexps
}