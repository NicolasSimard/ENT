/*ida=[a,(b+sqrt(D))/2], where a=N(ida), b is determined mod 2a and b^2=D mod 4a.*/
partial_Hecke(ida,l,s) = ida[1]^(-s)*imag(ida[2])^(2*l-s)*Eisen(l,ida[2]/ida[1],s)/2;

\r non-holo-Eisenstein.gp

myzetaK(D,s) =
{
    my(S=0,forms=reduced_forms(D));
    for(n=1,#forms,
        S += nonholoG(0,(forms[n][2]+sqrt(D))/2/forms[n][1],s);
    );
    (2/sqrt(abs(D)))^s*S/wD(D);
}
