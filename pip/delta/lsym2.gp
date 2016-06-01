/* Computing the Petersson norm of the delta function using the special values
of the symetric L-function attached to Delta. See Cohen Thm 2.1. Works. */

read("delta_LSym2f"); \\ Defines the symmetric square L-function (called L)

print("Petersson norm of Delta: ",L(12)*factorial(11)*2/Pi/(4*Pi)^12);