lsym2level1data(qexp, k, verb) = {
    my(a, Ldata);
    /*The formula for the coefficients of the L-function is taken from
    Cohen's article "Haberland's formula and numerical computation of the
    Petersson scalar product", Theorem 2.1.*/
    a = (n -> vector(n,i,sumdiv(i, d, (-1)^bigomega(d)*d^(k - 1)*qexp[i/d]^2)));

    Ldata = [a, 1, [0, 1, 2 - k], 2*k - 1 , 1, 1];
    
    if(verb, printf("Estimated error on symmetric square L-function: %.2f%%", (1.0-abs(lfuncheckfeq(lfuncreate(Ldata)))/default(realbitprecision))*100));
    
    return(Ldata);
}
addhelp(lsym2level1data, "lsym2level1data(qexp, k, {verb = 0}): Returns the L-function data attached to the modular form with q-expansion qexp (a vector), weigth k and level 1. If verb = 1, print the error in the functional equation (should be small!).");

sub(f, k, i) = {
    intnum(x = -1/2, 1/2,
        intnum(y = (1 - x^2)^(1/2),[[1],4 * Pi],
            norm(f((x + y * I)/(i * (x + y * I) + 1))*(i * (x + I * y) + 1)^-k)*y^(k - 2)));
}

sub0(f, k) = {
    intnum(x = -1/2, 1/2,
        intnum(y = (1 - x^2)^(1/2),[[1],4 * Pi],
            norm(f(-1/(x + y*I))*(x + I * y)^-k)*y^(k - 2)));
}

peterssonnorm(X, k, p, verb) = {
    if(!isprime(p) && p != 1,
        error("The level has to be prime of 1."));
    
    if(p == 1 && type(X) == "t_VEC", \\ Level one q-expansion, we use the Symmetric square L-function
        return((Pi / 2*(4 * Pi)^k/(k - 1)!)^-1*lfun(lfuncreate(lsym2level1data(X, k, verb)), k)));
    
    my(f);
    
    if(type(X) == "t_VEC", f = (z) -> sum(n = 1, #X, X[n]*exp(2 * Pi * I * n * z)), f = X);
    
    if(p != 1,
        if(verb, print(0));
        sub0(f, k))
    + sum(i = 0, p - 1, 
        if(verb, print((i + 1)/(p + 1)));
        sub(f, k, i));
}
addhelp(peterssonnorm, "peterssonnorm(X, k, p, {verb = 0}): compute the Petersson norm of the modular form given by the data X of weight k and level p, where p is prime or 1. Here X is either a function f from C to C (i.e. a closure (x) -> f(x)) or a q-expansion (i.e. a vector). If verb = 1, print the progress of the computations. In general, the computations are done using the definition as a double integral, but when the level is one and X is a q-expansion, it is assumed that the modular form is a newform and the symmetric square L-function is used (much more efficient).");

/*
A few examples:

gp > \p 50
gp > \ps 500
gp > \r pip/peterssonnorm.gp
gp > delta(z) = eta(z, 1)^24
gp > peterssonnorm((z) -> delta(z), 12, 1) \\ Using the def
%1 = 1.0353620568043209223478168122251645932249079609504 E-6
gp > peterssonnorm(Vec(delta(q)), 12, 1) \\ Using the symmetric square L-function
%2 = 1.0353620568043209223478168122251645932249079609373 E-6
gp > delta5(z) = (eta(z,1)*eta(5*z,1))^4
gp > peterssonnorm(x -> delta5(x), 4, 5)
%3 = 0.00087080010497868327553690871430136963039092614131510
gp > delta7(z) = (eta(z,1)*eta(7*z,1))^3 \\ delta7 has a non-trivial nebentype
gp > peterssonnorm(x->delta7(x), 3, 7)
%4 = 0.0052288338546578073366557434887148453071348136738030
gp > \\ An example where we start directly from the q-expansion of delta7
gp > qexp(pr) = Vec(q*prod(n = 1, pr, ((1 - q^n + O(q^pr)) * (1 - q^(7*n) + O(q^pr)))^3));
gp > delta7(z, pr) = my(a = qexp(pr)); sum(n = 1, #a, a[n]*exp(2*Pi*I*n*z));
gp > peterssonnorm((x) -> delta7(x, 10), 3, 7)
%5 = 0.0048706105843147222825918367692621948291339726670360
gp > peterssonnorm((x) -> delta7(x, 20), 3, 7)
%6 = 0.0056496863402561967596385041188041209738646270942585
*/