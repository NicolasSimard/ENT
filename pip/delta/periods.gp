/***Computing the periods of the delta modular form.                    ***/

/*Setting t_0 = 1 in the formula on page 13 of Cohen's paper.

The function uses suminf, since the terms decrease quickly. This means that the
summation is computed until the relative error of three consecutive terms is
less than the precision.
*/
delta_r(j) = {
    local(c);
    c = 2*Pi;
    I^(j+1)*suminf(n = 1, ramanujantau(n)*(incgam(j+1,c*n)/(c*n)^(j+1)\
                                    +incgam(12-j-1,c*n)/(c*n)^(12-j-1)))
}
