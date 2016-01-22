/*This script is meant to compute the q-expansions of all the theta series
attached to a quadratic field, up to a precision, with parameter k.
*/

read("../quadratic.gp");

\\ Returns the q-expansion of the theta function attached to the quadratic
\\ form with parameter k up to O(q^nbr_coeff).

theta_qexp(f,k,nbr_coeff) = {
    local(D, Bound, qexp, n, alpha, beta);

    if(k < 0, error("k has to be >= 0."));

    D = f[2]^2-4*f[1]*f[3];

    \\ Initialise the q-expansion with zeroes
    qexp = vector(nbr_coeff-1);

    \\ The generators of the ideal coresponding to f
    alpha = f[1];
    beta = (-f[2]+sqrt(D))/2;

    Bound = floor(sqrt(2*(nbr_coeff-1))/sqrt(f[1]+f[3]-sqrt((f[1]-f[3])^2+f[2]^2)));
    for(x = -Bound, Bound,
        for(y = -Bound, Bound,
            n = f[1]*x^2 + f[2]*x*y + f[3]*y^2;
            if(n < nbr_coeff && n != 0,
                qexp[n] += (alpha*x - beta*y)^(2*k);
            );
        );
    );
    \\ Note that qexp[n] = c_n, so there is no constant term
    qexp
}

\\ Returns the q-expansion of all the theta function attached to the quadratic
\\ field of discriminant D with parameter k up to O(q^nbr_coeff).

theta_qexps(D, k, nbr_coeff) = {
    local(forms);

    forms = reduced_forms(D);
    vector(length(forms),i,theta_qexp(forms[i],k,nbr_coeff))
}

eval_qexp(qexp,tau) = suminf(n=1,qexp[n]*exp(2*Pi*I*n*tau));

\\ This function computes the q-expansion at every evaluation, so it is highly
\\ inneficient if the function has to be evaluated many times (e.g. if it is
\\ integrated). Instead, one should compute the q-expansion once and then use
\\ eval_qexp... With 100 coefficients, one is usially able to evaluate up to
\\ 500 decimals at least.

theta_function(f,k,tau) = {
    local(qexp);
    eval_qexp(theta_qexp(f,k,100),tau)
}


