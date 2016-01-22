/* This script computes the individual entries a matrix where the (i,j)
entry is theta_i(tau_j)
*/
nbr_coeff = 100;

read("theta_qexp.gp")
read("../quadratic.gp");

\\ (i,j)th coefficient of the matrix consisting of [theta_i(tau_j)]_ij
a(D,k,i,j) = {
    local(exps,taus);
    exps = theta_qexps(D,k,nbr_coeff);
    taus = reduced_roots(D);
    suminf(n=1,exps[i][n]*exp(2*Pi*I*n*taus[j]))
}


eval_mat(D,k) = {
    local(exps,taus);
    exps = theta_qexps(D,k,nbr_coeff);
    taus = reduced_roots(D);
    matrix(length(taus),length(taus),i,j,suminf(n=1,exps[i][n]*exp(2*Pi*I*n*taus[j])))
}