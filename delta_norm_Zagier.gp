/*** Also trying using a formula for the special value of an        ***/
/*** L-function at 12 (see Zagier - Modular Forms).                 ***/

\\ This script defines the function tau(n).
read("delta_coeff.gp");

/* This function seems to work, but it converges very slowly and is long to
evaluate. N = 800 gives 9 digits and takes a few minutes to evaluate.*/

PeterssonNormZagier(N) = Pi*factorial(11)/(3*(4*Pi)^12)*sum(n=1,N,tau(n^2)/n^12)