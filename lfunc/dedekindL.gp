/* This script defines the Dedekind zeta function of a number field.

The input is simply an irreducible polynomial fpol defining the number field.
Then the weight, the sign, the gamma factor, the exponential factor and the
location of the pole are known. The residue of the completed L-function is
automatically computed.

Note: the residues at the poles are determined automatically. This gives an
efficient way of computing the class number, for example. However, if you only
want to evaluate the zeta function, it is better to use the built-in zetak()
function.

The function initDedekindL(fpol,verbose) resturns a vector [L,fullgamma,Lresidues],
where L*fullgamma is the completed L-function and Lresidues are the residues
of the COMPLETED L-function.
*/

initDedekindL(fpol = 1, verbose = 1) = {
    local(nf, r1, r2);
    nf        = nfinit(fpol);        \\ initialize the number field F/Q
    r1        = nf.r1;
    r2        = nf.r2;

    read("../computel");                     \\ read the ComputeL package

    conductor = abs(nf.disc);                \\ exponential factor
    weight    = 1;                           \\ L(s)=sgn*L(weight-s)
    sgn       = 1;                           \\ sign in the functional equation
    gammaV    =concat(vector(r1+r2,X,0),vector(r2,X,1));
    Lpoles    = [1];                         \\ pole at s=1, Lresidues=automatic
    \\This is the residue given by the class number formula (but times -1...)
    \\Lresidues = [-fullgamma(1)*2^bnf.r1*(2*Pi)^bnf.r2*bnf.reg*bnf.clgp.no/bnf.tu[1]/sqrt(abs(bnf.disc))];
    dzk       = dirzetak(nf,cflength());     \\ coefficients a(k) in L(s)

    initLdata("dzk[k]");         \\ initialize L-series

    \\Verify the functional equation. If verbose = 1, also print a message.
    \\This is also where the residue (whic is set to automatic) is computed.
    if(verbose,\
        print("Defining the zeta function of the number field defined by ",fpol);\
        print("Error in func. eq.  = ",errprint(checkfeq()));,checkfeq();\
    );
    [L,fullgamma,Lresidues]
};


