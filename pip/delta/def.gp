/***Computing the Petersson norm of delta directly from the definition.   ***/
/***This code gives the right answer!                                     ***/

\\delta(z) = eta(z,1)^24; \\The 1 indicates that we want the "true" eta (with q^(1/24))

print("Computing the norm using the definition.");
N = intnum(x = -1/2, 1/2, intnum(y = (1-x^2)^(1/2),[[1],4*Pi], norm(delta(x+y*I))*y^10));
print("The norm is: ", N);