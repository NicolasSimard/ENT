/* Computing the theta determinant attached to a quadratic field using the
definition of the Petersson inner product.
*/

\\ Define the parameters. Can be owerwritten.
p = 3;
Darmonk = 1;
nbr_coeff = 50;

print("Computing q-expansion.");
read("theta_qexp.gp");
qexps = theta_qexps(-p, Darmonk, nbr_coeff);

\\ The weight of the theta series
k = 2*Darmonk + 1;

\\ The class number of the quadratic field of discriminant -p
h = length(qexps);

\\ Define infinity
oo = [1];

f(i, z) = suminf(n=1, qexps[i][n]*exp(2*Pi*I*n*z));

sub(r, i, j) = intnum(x = -1/2, 1/2,intnum(y = (1-x^2)^(1/2),[oo,4*Pi],\
               f(i,(x+y*I)/(r*(x+y*I)+1))*conj(f(j,(x+y*I)/(r*(x+y*I)+1)))*\
               norm((r*(x+I*y)+1)^-k)*y^(k-2)));

sub0(i, j) = intnum(x = -1/2, 1/2,intnum(y = (1-x^2)^(1/2),[oo,4*Pi],\
             f(i, -1/(x+y*I))*conj(f(j, -1/(x+y*I)))*\
             norm((x+I*y)^-k)*y^(k-2)));

Petersson_inner(i, j) = (sub0(i, j) + sum(r = 0, p-1, sub(r, i, j)))/(p+1);

print("Computing the matrix for (D=",-p,", k=",Darmonk,"): ");

M = matrix(h,h);
{for(i=1,h,
    for(j=1,i,
        M[i,j] = Petersson_inner(i,j);
        M[j,i] = conj(M[i,j]);
        print(round((i*(i-1)/2+j)/(h*(h+1)/2)*100.),"%");
    );
);}

print("(D=",-p,", k=",Darmonk,"): ",matdet(M));