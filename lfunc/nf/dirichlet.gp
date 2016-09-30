/* This script defines the Dirichlet L-function of a PRIMITIVE Dirichlet char.

The input is a primitive Dirichlet character chi and its modulus (=conductor since
the chi is primitive).

The function initDedekindL(chi,modulus,verbose) returns a vector [L,fullgamma]
where L*fullgamma is the completed L-function.
*/

triv(n) = 1;

initDirichletL(chi = triv, modulus = 1, verbose = 1) = {
    local(a);

    read("../computel");

    conductor   = modulus;
    gammaV      = if(chi(-1) == 1,[0],[1]);  \\[0] if chi is even, [1] otherwise.
    weight      = 1;
    gausssum    = sum(k=0,modulus-1,chi(k)*exp(2*Pi*I*k/modulus));
    sgn         = gausssum/I^gammaV[1]/sqrt(modulus);

    Lpoles      = if(modulus == 1, [1],[]);
    Lresidues   = if(modulus == 1, [-1],[]);

    a = vector(cflength(),k,chi(k));
    initLdata("a[k]",,"conj(a[k])");

    \\Verify the functional equation. If verbose = 1, also print a message.
    if(verbose,\
        print("Defining the Dirichlet L-function of a character of modulus ",modulus);\
        print("Error in func. eq.  = ",errprint(checkfeq()));\
        print("The sign is         = ",sgn),\
        checkfeq(););
    [L,fullgamma]
}

