/* Computing the theta determinant attached to a quadratic field using the
definition of the Petersson inner product.
*/

\rtheta_qexp.gp;

eval_form(coeff, z) = suminf(n=1, coeff[n]*exp(2*Pi*I*n*z));

sub(r, f1, f2, k) = intnum(x = -1/2, 1/2,intnum(y = (1-x^2)^(1/2),[[1],4*Pi],\
               eval_form(f1,(x+y*I)/(r*(x+y*I)+1))*\
               conj(eval_form(f2,(x+y*I)/(r*(x+y*I)+1)))*\
               norm((r*(x+I*y)+1)^-k)*y^(k-2)));

sub0(f1, f2, k) = intnum(x = -1/2, 1/2,intnum(y = (1-x^2)^(1/2),[[1],4*Pi],\
             eval_form(f1, -1/(x+y*I))*conj(eval_form(f2, -1/(x+y*I)))*\
             norm((x+I*y)^-k)*y^(k-2)));

Petersson_inner(f1, f2, k, p) = (sub0(f1, f2, k) + sum(r = 0, p-1, sub(r, f1, f2, k)))/(p+1);

theta_det(p, Darmonk, verbose = 0) = {
    if(verbose, print("Computing q-expansion."));
    nbr_coeff = 100;
    qexps = theta_qexps(-p, Darmonk, nbr_coeff);

    \\ The weight of the theta series
    k = 2*Darmonk + 1;

    \\ The class number of the quadratic field of discriminant -p
    h = length(qexps);
    M = matrix(h,h);
    if(verbose, print("Computing the matrix."));
    for(i=1,h,
        for(j=1,i,
            M[i,j] = Petersson_inner(qexps[i],qexps[j], k, p);
            M[j,i] = conj(M[i,j]);
            if(verbose, print(round((i*(i-1)/2+j)/(h*(h+1)/2)*100.),"%"));
        );
    );
    matdet(M)
};

