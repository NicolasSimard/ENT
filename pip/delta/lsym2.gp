/*Computing the Petersson norm of delta using the special values of its symmetric square
L-function.

Tested up to the number of decimals in Cohen's paper on Haberland's formulas.

Here the Petersson inner product is defined as \int\int_{\SL_2(Z)\H}...dxdy/y^2,
i.e. no normalization.
*/

k = 12; \\ Weight of delta
N = 1;

\\ Initialize data for the L-function L(sym^2(delta),s)
a = (n -> vector(n,i,sumdiv(i, d, (-1)^bigomega(d)*d^(k-1)*ramanujantau(i/d)^2)));
astar = 1;
Vga = [0,1,2-k];
w = 2*k-1;
cond = N^2;
eps = 1;

Ldata = [a,astar,Vga,w,cond,eps];

printf("Estimated error on symmetric square L-function: %.2f%%", (1.0-abs(lfuncheckfeq(lfuncreate(Ldata)))/default(realbitprecision))*100);

print("The Petersson norm of delta is (stored as delta_norm): ");
delta_norm = (Pi/2*(4*Pi)^k/(k-1)!/N)^-1*lfun(lfuncreate(Ldata),k)