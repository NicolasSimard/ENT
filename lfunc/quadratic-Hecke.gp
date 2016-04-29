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

/* Hecke L-function of a Hecke character char of conductor 1 and infinity type t, i.e.
char:I_K--->C^\times is a group homomorphism from the group of fractional ideals
to C^\times such that psi((x))=x^t for all x in K. Note that when disc(K) < -4,
such a character is well-defined if and only if t=2l is even.

char is given as a list of pairs [ida,psi(ida)], one for each representative ida of ClK.
Note that ida has to be given as a Z-basis [a,(b+sqrt(D))/2], where a=N(ida), b is
determined mod 2a and b62=D mod 4a.
*/
quadHeckeL(char,t,s) =
{
    my(S,sqrtD);
    \\char[n] = [ida,psi(ida)]
    S = sum(n=1,#char,char[n][2]/char[n][1][1]^t*nonholoG(t,char[n][1][2]/char[n][1][1],s-t));
    sqrtD = 2*imag(char[1][1][2]);
    (2/sqrtD)^(s-t)*S/wD(-round(sqrtD^2));
}
addhelp(quadHeckeL,"quadHeckeL(char,t,s): char=[[ida1,chi(ida1)],...,[idah,psi(idah)]],t=type.")

trivHecke(D) =
{
    my(forms=reduced_forms(D));
    vector(#forms,n,[ida(forms[n]),1]);
}
