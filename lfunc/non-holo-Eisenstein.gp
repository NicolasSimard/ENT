MiyakeA(k,s) = 2^(k+1)*I^-k*Pi^(s+k)*if(type(s) == "t_INT" && s+k <= 0,0,gamma(s+k)^-1);

MiyakeB(k,s) = 2^(1-k)*I^-k*Pi^s*if(type(s) == "t_INT" && s <= 0,0,gamma(s)^-1);

MiyakeC(k,s) = 2*zeta(2*s+k);

GProd(k,s) =
{
    if(type(s) == "t_INT",
        if(s + k   <= 0, return(0));
        if(s + k/2 <= 0, return((-1)^(k/2)*gamma(s)/gamma(n+k)/gamma(n+k/2+1)));
        if(s       <= 0, return(0));
    );
    return(gamma(s+k/2)/gamma(s+k)/gamma(s));
}

zetastar(s) = Pi^(-s/2)*gamma(s/2)*zeta(s);

MiyakeD(k,s) = 2*Pi^(s+k/2)*I^-k*if(s == 1-k/2,4*sqrt(Pi)*(-1)^(k/2-1)/k,GProd(k,s)*zetastar(2*s+k-1));

Miyakean(n,k,s) = sigma(n,k+2*s-1);

MiyakeSum1(k,z,s) = (4*Pi*imag(z))^s*suminf(n=1,Miyakean(n,k,s)*exp(2*Pi*I*n*z)*hyperu(s,k+2*s,4*Pi*imag(z)*n));

MiyakeSum2(k,z,s) = (4*Pi*imag(z))^(k+s)*suminf(n=1,Miyakean(n,k,s)*exp(2*Pi*I*n*-conj(z))*hyperu(k+s,k+2*s,4*Pi*imag(z)*n));

MiyakeE(k,z,s) =
{
    MiyakeC(k,s) + MiyakeD(k,s)*imag(z)^(1-k-2*s)
    + MiyakeA(k,s)*imag(z)^-s*MiyakeSum1(k,z,s)
    + MiyakeB(k,s)*imag(z)^(-s-k)*MiyakeSum2(k,z,s);
}
addhelp(MiyakeE,"MiyakeE(k,z,s): Non-holomorphic Eisenstein series defined as in Miyake's book"\
" book on modular forms for the trivial characters, i.e."\
" Miyake(k,z,s) = sum_{m,n}(mz+n)^(-k)|mz+n|^(-2s).")

nonholoG(k,z,s) = imag(z)^s*MiyakeE(k,z,s);
addhelp(nonholoG,"nonholoG(k,z,s): Non-holomorphic Eisenstein series defined as in Miyake's book"\
" book on modular forms for the trivial characters, BUT multiplied by y^s, i.e."\
" nonholoG(k,z,s) = Im(z)^s*sum_{m,n}(mz+n)^(-k)|mz+n|^(-2s).")
