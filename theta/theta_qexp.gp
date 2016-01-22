/*This script is meant to compute the q-expansions of all the theta series
attached to a quadratic field, up to a precision, with parameter k.
*/

read("../quadratic.gp");

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
