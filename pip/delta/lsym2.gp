/*Computing the Petersson norm of delta using the special values of its symmetric square
L-function.

Tested up to the number of decimals in Cohen's paper on Haberland's formulas.

Here the Petersson inner product is defined as \int\int_{\SL_2(Z)\H}...dxdy/y^2,
i.e. no normalization.
*/

w=12; \\ Weight of delta

\\ Initialize data for the L-function L(sym^2(delta),s)
a = (n -> vector(n,i,sumdiv(i, d, (-1)^bigomega(d)*d^(w-1)*ramanujantau(i/d)^2)));
astar = 1;
Vga = [0,1,2-w];
k = 2*w-1;
N = 1;
eps = 1;

\\ print("Check functionnal equation (should be very negative...):", lfuncheckfeq(lfuncreate([a,astar,Vga,k,N,eps])));
print("The Petersson norm of delta is (stored as delta_norm): ");
delta_norm=2/Pi*factorial(w-1)/(4*Pi)^w*lfun(lfuncreate([a,astar,Vga,k,N,eps]),w)
