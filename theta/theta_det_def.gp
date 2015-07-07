/*** Computing the discriminant of the Petersson inner product on the     ***/
/*** theta series attached to the field of discriminant -3                ***/

\\ Define the parameters
p = 3;
Darmonk = 3;
nbr_coeff = 50;

print("Computing q-expansion");
\\ Compute the q-expansion and store it in a file as a list called qexps
output = concat("qexps/",concat(Str(p),concat("_",concat(Str(Darmonk),concat("_",concat(Str(nbr_coeff),".gp"))))));
commande = concat("python theta_qexp.py ",concat(Str(-p),concat(" ",concat(Str(Darmonk),concat(" ",concat(Str(nbr_coeff)," pari"))))));
system(concat(commande,concat(" > ",output)));

\\ Read the q-expansion. This defines qExps
read(output);

\\ The keight of the theta series
k = 2*Darmonk + 1;

\\ The class number of the quadratic field of discriminant -p
h = length(qExps);

\\ Define infinity
oo = [1];

f(i, z) = suminf(n=1, qExps[i][n]*exp(2*Pi*I*n*z));

sub(r, i, j) = intnum(x = -1/2, 1/2,intnum(y = (1-x^2)^(1/2),[oo,4*Pi],\
               f(i,(x+y*I)/(r*(x+y*I)+1))*conj(f(j,(x+y*I)/(r*(x+y*I)+1)))*\
               norm((r*(x+I*y)+1)^-k)*y^(k-2)));

sub0(i, j) = intnum(x = -1/2, 1/2,intnum(y = (1-x^2)^(1/2),[oo,4*Pi],\
             f(i, -1/(x+y*I))*conj(f(j, -1/(x+y*I)))*\
             norm((x+I*y)^-k)*y^(k-2)));

Peter_inner(i, j) = {
    if(i > j,
        conj(Petersson_inner(j,i)),
        (sub0(i, j) + sum(r = 0, p-1, sub(r, i, j)))/(p+1)
    );
};

M = matrix(h,h);

print("Computing the matrix...");
{for(i=1,h,
    for(j=1,i,
        M[i,j] = Peter_inner(i,j);
        M[j,i] = conj(M[i,j]);
        print(round((i*(i-1)/2+j)/(h*(h+1)/2)*100.),"%");
    );
);}

print("(D=",-p,", k=",Darmonk,"): ",matdet(M));