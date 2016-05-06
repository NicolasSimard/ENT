/* This script uses he formula for the Petersson norm of theta series in terms
of linear combinations of derivatives of the Eisenstein series E2.*/

E2num(z) = (8*Pi*imag(z))^-1-1/24+suminf(n=1,sigma(n)*exp(2*Pi*I*n*z));

E4num(z) = 1/240+suminf(n=1,sigma(n,3)*exp(2*Pi*I*n*z));

E6num(z) = -1/504+suminf(n=1,sigma(n,5)*exp(2*Pi*I*n*z));

dnE2(n) =
{
    my(v,d);
    v=['E2,'E4,'E6];
    d=[5/6*'E4-2*'E2^2,7/10*'E6-8*'E2*'E4,400*'E4^2-12*'E2*'E6];
    z->subst(subst(subst(diffop('E2,v,d,n),'E2,E2num(z)),'E4,E4num(z)),'E6,E6num(z));
}

dn(P,n) =
{
    my(v,d);
    v=['E2,'E4,'E6];
    d=[5/6*'E4-2*'E2^2,7/10*'E6-8*'E2*'E4,400*'E4^2-12*'E2*'E6];
    diffop(P,v,d,n);
}

dnnum(P,n) =
{
    my(v,d);
    v=['E2,'E4,'E6];
    d=[5/6*'E4-2*'E2^2,7/10*'E6-8*'E2*'E4,400*'E4^2-12*'E2*'E6];
    z->subst(subst(subst(diffop(P,v,d,n),'E2,E2num(z)),'E4,E4num(z)),'E6,E6num(z));
}