/*** Computing the Petersson norm of the delta function using the special ***/
/*** values of the symetric L-function attached to Delta. See Cohen p.4.  ***/

\\ This script defines the function tau(n).
read("delta_coeff.gp");

\\The weight, to simplify
k = 12;

\\The fourrier coefficients of L(sym2f,s)
A(n) = sumdiv(n, d, (-1)^bigomega(d)*d^(k-1)*tau(n/d)^2);

\\A constant needed in the formula
C = 2*Pi^(3/2);

\\The harmonic series
H(n) = sum(k=1,n,1./k);

\\a gamma factor needed in the formula
g(s) = C^(-s)*gamma(s)*gamma((s-k)/2+1);

\\The components of the function

F1(s,x) = sum(m=1, (k-2)/2, (-1)^(k/2-m-1)*factorial(2*m-1)/factorial(k/2-m-1)*(C*x)^(-2*m)/(s-2*m));

F2(s,x) = 2^(k-1)*suminf(m=0, (-1)^(k/2-m-1)*factorial(m+k/2)/factorial(2*m+1)/factorial(2*m+k)/(s+2*m+1)*(2*C*x)^(2*m+1));

F3(s,x) = suminf(m=0, (-1)^(k/2-m-1)*(C*x)^(2*m)/factorial(2*m)/factorial(m+k/2+1)/(2*m+s)*\
          (2*H(2*m)+H(m+k/2-1)-3*Euler-2*log(C*x)+2./(2*m+s)));

\\The main function
F(s,x) = g(s)-x^s*(2*F1(s,x)+Pi^(1/2)*F2(s,x)+F3(s,x));

print("The computed inner product is: ");
Norm = 2^(1-k)*Pi^(k/2-1)*suminf(n=1, A(n)/n^k*(F(k,n)-n*F(k-1,n)));
print(Norm);