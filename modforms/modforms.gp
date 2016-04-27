/* This script defines a few modular forms. These modular forms are represented
as functions sending integers n to a_n(f), the nth Fourier  coefficient of f.

Some modular forms can also be represented as Fourrier series directly. This is
often faster than computing every Fourrier coefficient one after the other. For
example, delta_qexp(1000) is much faster than Ser(vector(1000,n,delta(n))).*/

\r operators.gp

clos2qexp(f,pr:small) = Ser(vector(pr,n,f(n-1)),'q);
addhelp(clos2qexp,"cols2qexp(f,pr): Return the q-expansion attached to the closure representation of f up to precision pr.");

qexp2clos(f) = n -> polcoeff(f,n,'q);
addhelp(qexp2clos,"qexp2clos(f): Returns the closure attached to the q-expansion of f.");

/* To define the Eisenstein series, we follow Zagier's convention:

G_k(s) = 1/2*sum_{(m,n)\neq(0,0)}(mz+n)^-k

E_k(s)=G_k(s)/zeta(k) (leading coefficient is 1)

GG_k(s) = (k-1)!/(2*Pi*I)^kG_k(s) (all non-constant Fourrier coeff are int)
*/

GG(k:small) = n -> if(n == 0, -bernfrac(k)/2/k, sigma(n,k-1));
addhelp(GG,"GG(k): Returns Eisenstein series GG of weight k. Modular forms are\nrepresented as closures, sending n to the nth Fourrier coefficient.");

GG_qexp(k:small,pr:small) = -bernfrac(k)/2/k+Ser(concat([0],vector(pr-1,n,sigma(n,k-1))),'q);
addhelp(GG_qexp,"GG_qexp(k,pr): Returns the q-expansion of the Eisenstein series\nGG of weight k.");

E(k:small) = n -> if(n == 0, 1, -2*k/bernfrac(k)*sigma(n,k-1));
addhelp(E,"E(k): Returns Eisenstein series E of weight k. Modular forms are\nrepresented as closures, sending n to their nth Fourrier coefficient.");

E_qexp(k:small,pr:small) = -2*k/bernfrac(k)*GG_qexp(k,pr);
addhelp(E_qexp,"E_qexp(k,pr): Returns the q-expansion of the Eisenstein series\nE of weight k.");

eta3_qexp(pr:small) =  my(d); Ser(vector(pr,n,if(issquare(1+8*(n-1),&d),(-1)^((-1+d)/2)*d,0)),'q);

theta_qexp(pr:small) = Ser(concat([0],vector(pr-1,n,2*issquare(n))),'q)+1;
addhelp(theta_qexp,"theta_qexp(k): Returns the q-expansion of the theta series (=sum_{n\in\Z} q^(n^2)) of weight 1/2 up to precision pr.");

mftheta(n) = if(n==0, 1, if(issquare(n),2,0));
addhelp(mftheta,"mftheta(n): nth Fourrier coefficient of the theta series (=sum_{n\in\Z} q^(n^2)) of weight 1/2");

theta1_qexp(pr:small) = my(d); Ser(concat([0],vector(pr-1,n,2*issquare(n,&d)*(-1)^d)),'q)+1;
addhelp(theta1_qexp,"theta1_qexp(k): Returns the q-expansion of the theta series (=sum_{n\in\Z} (-1)^n*q^(n^2)) of weight 1/2 up to precision pr.");

mftheta1(n) = my(d); if(n==0, 1, if(issquare(n,&d),2*(-1)^d,0));
addhelp(mftheta1,"mftheta1(n): nth Fourrier coefficient of the theta1 series (=sum_{n\in\Z} (-1)^n*q^(n^2)) of weight 1/2");

delta(n) = if(version[2] >= 8, print("ram"); ramanujantau(n), polcoeff(delta_qexp(n+1),n,'q));
addhelp(delta,"delta: Returns the modular form delta. Modular forms are represented\nas functions sending n to their nth Fourrier coefficient.");

delta_qexp(pr:small) = 'q*sqr(sqr(sqr(eta3_qexp(pr-1))));
addhelp(delta_qexp,"delta_qexp(k): Returns the q-expansion of the delta series up to precision pr.");

j_qexp(pr:small) = E_qexp(4,pr+2)^3/delta_qexp(pr+2);
addhelp(j_qexp,"j_qexp(pr): Returns the q-expansion of the j-invariant up to precision pr.");

jpol(f) =
{
    my(j,js,M,B,k=-valuation(f,'q));
    if(k <=  0, return(f));
    j=j_qexp(k);
    js=vector(k+1,n,j^(k+1-n));
    M=matrix(k+1,k+1,n,m,polcoeff(js[m],n-k-1,'q));
    B=vector(k+1,n,polcoeff(f,n-k-1,'q));
    return(Pol(matsolve(M,B~),'X));
}
addhelp(jpol,"jpol(f): Given a q-expansion f(q) with enough coefficients, returns a polynomial g(Y) such that g(j)=f, where j is the j-function.");

clos2qexp(f,pr:small) = Ser(vector(pr,n,f(n-1)),'q);
addhelp(clos2qexp,"cols2qexp(f,pr): Return the q-expansion attached to the closure representation of f up to precision pr.");

qexp2clos(f) = n -> polcoeff(f,n,'q);
addhelp(qexp2clos,"qexp2clos(f): Returns the closure attached to the q-expansion of f.");

fi(i:small, pr:small = 15) =
{
    if(i%4 != 0 && i%4 != 3, error("Wrong value for i: has to be cong to 0 or 3 mod 4."));
    my(fis = vector(i+1), tmp, j4); \\ fis[i] = f_{i-1}, so f_i = fis[i+1]
    if(i == 0, return(theta_qexp(pr)), fis[1] = theta_qexp(pr));
    tmp = (theta_qexp(pr)*d(V(4)(E_qexp(10,pr)))/2-10*d(theta_qexp(pr))*V(4)(E_qexp(10,pr)))/V(4)(delta_qexp(pr));
    tmp = (tmp+608*fis[1])/-20;
    if(i == 3, return(tmp), fis[4] = tmp);
    j4 = V(4)(j_qexp(pr));
    for(j = 4, i,
        if(j%4 == 0 || j%4 == 3,
            tmp = j4*fis[j-4+1];
            for(n = -j + 1, 0,
                if((-n)%4 == 0 || (-n)%4 == 3, 
                    tmp -= polcoeff(tmp,n,'q)*fis[-n+1]
                );
            );
            fis[j+1] = tmp;
        );
    );
    fis[i+1];
}

fi2(i:small, pr:small = 20) =
{
    if(i%4 != 0 && i%4 != 3, error("Wrong value for i: has to be cong to 0 or 3 mod 4."));
    my(g1, g1quot, g4, g4quot, j4, j4pow, f0, f3, fi = 0);
    if(i == 0, return(theta_qexp(pr)), f0 = theta_qexp(pr));
    f3 = (theta_qexp(pr)*d(V(4)(E_qexp(10,pr)))/2-10*d(theta_qexp(pr))*V(4)(E_qexp(10,pr)))/V(4)(delta_qexp(pr));
    f3 = (f3+608*f0)/-20;
    if(i == 3, return(f3));
    
    j4 = V(4)(j_qexp(pr));
    g1 = theta1_qexp(pr)*V(4)(E_qexp(4,pr))/V(4)(eta3_qexp(pr))^2/'q;
    g4 = (10*d(g1)*V(4)(E_qexp(10,pr))-3/2*g1*d(V(4)(E_qexp(10,pr))))/V(4)(delta_qexp(pr));
    g4 = (g4 + (10*j4-21344)*g1)/-20;

    j4pow = 1;
    g1quot = g1/j4;
    g4quot = g4/j4;

    for(n =0, i\4 - (i%4 == 0),
        fi += polcoeff(g4quot,i,'q)*j4pow*f0 + polcoeff(g1quot,i,'q)*j4pow*f3;
        j4pow *= j4;
        g1quot /= j4;
        g4quot /= j4;
    );
    fi += polcoeff(g4quot,i,'q)*j4pow*f0;
}
