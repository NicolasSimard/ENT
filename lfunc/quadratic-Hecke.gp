/*Terms of the Fourier expansion of the non-holomorphic Eisenstein series
G_2l(z,s-2l).*/

term1(l,z,s) = 2*imag(z)^(s-2*l)*zeta(2*s-2*l);

zetastar(s) = Pi^(-s)*gamma(s/2)*zeta(s);

term2(l,z,s) = 
{
    if(s <= 0, return(0));
    2*(-1)^l*imag(z)^(s-2*l)*Pi^(s-l)/gamma(s)*prod(k=1,l,s-l-k)*zetastar(2*s-2*l-1);
}

term3(l,z,s) =
{
    my(S);
    if(s <= 0, return(0));
    S = suminf(n=1,sigma(n,2*s-2*l-1)*exp(2*Pi*I*n*z)*hyperu(s-2*l,2*s-2*l,4*Pi*imag(z)*n));
    2*(-1)^l*(2*Pi)^(2*s-2*l)*imag(z)^(s-2*l)/gamma(s)*S;
}

term4(l,z,s) =
{
    my(S);
    if(s <= 2*l, return(0));
    S = suminf(n=1,sigma(n,2*s-2*l-1)*exp(-2*Pi*I*n*conj(z))*hyperu(s,2*s-2*l,4*Pi*imag(z)*n));
    2*(-1)^l*(2*Pi)^(2*s-2*l)*imag(z)^(s-2*l)/gamma(s-2*l)*S;
}

\\Note: Eisen(l,z,s)=G_2l(z,s-2l)
Eisen(l,z,s) = term1(l,z,s)+term2(l,z,s)+term3(l,z,s)+term4(l,z,s);

/*
ida=[a,(b+sqrt(D))/2], where a=N(ida), b is determined mod 2a and b^2=D mod 4a.
*/
partial_Hecke(ida,l,s) = ida[1]^(-s)*imag(ida[2])^(2*l-s)*Eisen(l,ida[2]/ida[1],s)/2;

myzetaK(D,s) =
{
    my(S=0,forms=reduced_forms(D));
    for(n=1,#forms,
        S += Eisen(0,tau(forms[n]),s);
    );
    (2/sqrt(-D))^s*S/2;
}
