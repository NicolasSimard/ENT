/*** Computing with the theta function attached to the ideal class      ***/
/*** correspondign to x^2+xy+6y^2 in the quadratic field of             ***/
/*** discriminant -23.                                                  ***/
/*** type \rtheta_norm.gp or read("theta_norm.gp") at pari prompt to run this     ***/

read("computel");                 \\ read the ComputeL package
                                  \\ and set the default values
default(realprecision,150);       \\ set working precision; used throughout

N = 23;
k = 1;

                            \\ initialize L-function parameters
conductor = sqrt(N)/Pi;     \\ exponential factor
gammaV    = [0,1];          \\ list of gamma-factors
weight    = 2*k+1;          \\ L(s)=sgn*L(weight-s)
sgn       = 1;              \\ sign in the functional equation (NOT SURE))
Lpoles    = [];

print("We need ",cflength()," coefficients.");

a = [0, 2, 0, 0, 8, 0, -14-8*(1 + sqrt(-23))/2, 0, -6-8*(1 + sqrt(-23))/2, 18, 0, 0, 10-8*(1 + sqrt(-23))/2, 0, 0, 0, 32, 0, 34-8*(1 + sqrt(-23))/2, 0, 0, 0, 0, -30-16*(1 + sqrt(-23))/2, -56-32*(1 + sqrt(-23))/2, 50, 66-8*(1 + sqrt(-23))/2, -44-32*(1 + sqrt(-23))/2, 0, 0, 0, 0, -24-32*(1 + sqrt(-23))/2, 0, 0, 0, 178-8*(1 + sqrt(-23))/2, 0, 0, 4-32*(1 + sqrt(-23))/2, 0, 0, 0, 0, 0, 0, 0, 0, 194-40*(1 + sqrt(-23))/2, 98, 0, 0, -134-72*(1 + sqrt(-23))/2, 0, -126-72*(1 + sqrt(-23))/2, 0, 0, 0, -110-72*(1 + sqrt(-23))/2, 84-32*(1 + sqrt(-23))/2, 0, 0, 210-8*(1 + sqrt(-23))/2, 0, 42-72*(1 + sqrt(-23))/2, 0, 0, 0, 0, 0, 0, 0, 82-104*(1 + sqrt(-23))/2, 0, 0, 0, 0, 0, 274-8*(1 + sqrt(-23))/2, 0, 0, 162, -14-72*(1 + sqrt(-23))/2, 0, 0, 0, 0, 196-32*(1 + sqrt(-23))/2, 0, 0, 0, 0, -120-64*(1 + sqrt(-23))/2, -236-128*(1 + sqrt(-23))/2, 34-72*(1 + sqrt(-23))/2, 0, 122-136*(1 + sqrt(-23))/2, 0, 0, 0];

initLdata("a[k]");          \\Has to be a[k]... a[n] doesn't work.

/*
period(n)   =(-2*Pi*I)^(-n-1)*factorial(n)*L(n+1);    \\The n-th period of Delta
Ind(m,n)    =if(Mod(m,2) != Mod(n,2),1,0);            \\Indicator function
S = sum(m=0,10,sum(n=0,10-m,(-1)^m*Ind(m,n)*binomial(10,m+n)*binomial(m+n,m)*period(m)*conj(period(n))));*/

print("Verifying functional equation. Error: ",errprint(checkfeq()));
print("Computing the Petterson norm of theta:");
print("1) Using Haberland's formula:   ",N/(3*(-2*I)^11));
