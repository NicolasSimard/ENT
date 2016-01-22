/* Computing the theta determinant attached to a quadratic field using the
definition of the Petersson inner product.
*/
nbr_coeff = 100;

read("theta_function.gp");

sub(r, f1, f2, k) = intnum(x = -1/2, 1/2,intnum(y = (1-x^2)^(1/2),[[1],4*Pi],\
               eval_qexp(f1,(x+y*I)/(r*(x+y*I)+1))*\
               conj(eval_qexp(f2,(x+y*I)/(r*(x+y*I)+1)))*\
               norm((r*(x+I*y)+1)^-k)*y^(k-2)));

sub0(f1, f2, k) = intnum(x = -1/2, 1/2,intnum(y = (1-x^2)^(1/2),[[1],4*Pi],\
             eval_qexp(f1, -1/(x+y*I))*conj(eval_qexp(f2, -1/(x+y*I)))*\
             norm((x+I*y)^-k)*y^(k-2)));

\\ Here k is the weight (not the parameter of the theta series...)

petersson_inner(f1, f2, k, p) = (sub0(f1, f2, k) + sum(r = 0, p-1, sub(r, f1, f2, k)))/(p+1);

\\ D is the discriminant of a quadratic field and k is the parameter of the
\\ theta series (so that the weight is 2*k+1).

petersson_mat(D, k, verbose = 0) = {
    local(h,M);

    if(isprime(-D) == 0, error("Not yet implemented for non prime discriminant"));

    qexps = theta_qexps(D, k, nbr_coeff);

    \\ The class number of the quadratic field of discriminant D
    h = length(qexps);
    M = matrix(h,h);
    if(verbose, print("Computing the matrix."));
    for(i=1,h,
        for(j=1,i,
            M[i,j] = petersson_inner(qexps[i],qexps[j], 2*k+1, -D);
            M[j,i] = conj(M[i,j]);
            if(verbose, print(round((i*(i-1)/2+j)/(h*(h+1)/2)*100.),"%"));
        );
    );
    M
};

