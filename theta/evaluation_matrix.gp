/* This script computes the individual entries a matrix where the (i,j)
entry is theta_i(tau_j)
*/
nbr_coeff = 100;

read("theta_function.gp")
read("../quadratic.gp");

\\ (i,j)th coefficient of the matrix consisting of [theta_i(tau_j)]_ij
a(D,dk,i,j) = {
    local(exps,taus);
    exps = theta_qexps(D,dk,nbr_coeff);
    taus = reduced_roots(D);
    suminf(n=1,exps[i][n]*exp(2*Pi*I*n*taus[j]))
}

\\Normalized coefficient. Has to be algebraic.
a_normal(D,dk,i,j) = {
    local(exps,taus);
    exps = theta_qexps(D,dk,nbr_coeff);
    taus = reduced_roots(D);
    suminf(n=1,exps[i][n]*exp(2*Pi*I*n*taus[j]))/CSperiod(D)^(2*dk+1)
}

eval_mat(D,dk) = {
    local(exps,taus);
    exps = theta_qexps(D,dk,nbr_coeff);
    taus = reduced_roots(D);
    matrix(length(taus),length(taus),i,j,suminf(n=1,exps[i][n]*exp(2*Pi*I*n*taus[j])))
}

\\Normalized matrix. Each entry hasto be algebraic.
eval_mat_normal(D,dk) = {
    local(exps,taus);
    exps = theta_qexps(D,dk,nbr_coeff);
    taus = reduced_roots(D);
    matrix(length(taus),length(taus),i,j,\
           suminf(n=1,exps[i][n]*exp(2*Pi*I*n*taus[j]))/CSperiod(D)^(2*dk+1))
}



