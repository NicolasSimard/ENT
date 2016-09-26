/*This script is meant to compute the q-expansions of all the theta series
attached to a quadratic field, up to a precision, with parameter dk.
*/

read("../quadratic.gp");

/* Returns the q-expansion of the theta function attached to the quadratic
 form with parameter dk up to O(q^nbr_coeff). The different formats are
 "list" (default) : return a list [c_0,c_1,c_2,...,c_(nbr_coeff-1)]
 "xs"             : return a list [l1,l2,...,l_(nbr_coeff-1)], where
                    ln=[x1,...,xk] s.t. N(xi)/N(f) = n.
 "series"         : return a q-expansion c_1*q+c_2*q^2+...
*/
theta_qexp(f,dk,nbr_coeff,format="list",exact=0) = {
    local(D, Bound, xs, list, n, alpha, beta);

    if(dk < 0, error("k has to be >= 0."));

    D = f[2]^2-4*f[1]*f[3];

    \\ Initialise the q-expansion with zeroes
    xs = vector(nbr_coeff-1,n,[]);

    \\ The generators of the ideal coresponding to f
    alpha = f[1];
    beta = if(exact,(-f[2]+rD)/2,(-f[2]+sqrt(D))/2);

    Bound = floor(sqrt(2*(nbr_coeff-1))/sqrt(f[1]+f[3]-sqrt((f[1]-f[3])^2+f[2]^2)));
    for(x = -Bound, Bound,
        for(y = -Bound, Bound,
            n = f[1]*x^2 + f[2]*x*y + f[3]*y^2;
            if(n < nbr_coeff && n != 0,
                xs[n] = concat(xs[n],[alpha*x - beta*y])
                \\qexp[n] += (alpha*x - beta*y)^(2*dk);
            );
        );
    );
    if(dk == 0,xs=concat([[0]],xs)); \\ xs now has length (nbr_coeff-1) + 1
    if(format == "xs", return(xs));
    list=vector(nbr_coeff,n,sum(i=1,length(xs[n]),xs[n][i]^(2*dk)));
    if(exact,list=lift(Mod(list,rD^2-D)));
    if(format == "list", return(list));
    if(format == "series",return(Ser(list,q)));
    error("Unknown format: ",format);
}

\\ Returns the q-expansion of all the theta function attached to the quadratic
\\ field of discriminant D with parameter dk up to O(q^nbr_coeff).

theta_qexps(D, dk, nbr_coeff) = {
    local(forms);

    forms = reduced_forms(D);
    vector(length(forms),i,theta_qexp(forms[i],dk,nbr_coeff))
}

eval_qexp(qexp,tau) = sum(n=0,length(qexp)-1,qexp[n+1]*exp(2*Pi*I*n*tau));

theta_function(f,dk,tau,prec) = {
    local(qexp);
    eval_qexp(theta_qexp(f,dk,prec),tau);
}
