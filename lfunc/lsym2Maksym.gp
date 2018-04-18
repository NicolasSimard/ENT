k = 12; \\ Weight of delta
N = 1;

\\ Initialize data for the L-function L(sym^2(delta),s)
\\ Formula taken from "Cohen - Haberland formulas..." Thm.2.1
a = (n -> vector(n,i,sumdiv(i, d, (-1)^bigomega(d)*d^(k-1)*ramanujantau(i/d)^2)))

L(A, delta, k, N, sym2f) = {
    my(x = (k * N * delta)^2);

    coeffs = if(type(coeffs) == "t_CLOSURE",
        print("computing ", x," coefficients");
        sym2f(x)
    ,
        sym2f);

    print("computing the sum");
    sum(n=1, x, coeffs[n]/n^((k+1)/2)*exp(-delta*n/x)^A)
}

/*
TO RUN THIS SCRIPT, GO TO THE ENT DIRECTORY AND LOAD IT:
gp > read("lfunc/lsym2Maksym.gp")
---------------------
Example 1: use the function a directly and it will compute the right amount of terms
gp > L(2, 10, 12, 1, a)
computing 14400 coefficients
computing the sum
%27 = 66939818812.87144165805278406755752943799536675544152334824011662970866380857037887985941042181415802
----------------------
Example 2: compute and store the coefficients (make sure to compute enough!)
gp > coeffs = a(14400)
%28 = [1, -1472[+++]
gp > L(3, 10, 12, 1, coeffs)
computing the sum
%29 = 3054557.76135892997684177721411744[...]
*/