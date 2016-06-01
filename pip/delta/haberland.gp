\\Define the function delta_r(n), which computes the periods of Delta.
read("delta_periods.gp");

Ind(m,n) = if(Mod(m,2) != Mod(n,2),1,0);

print("Computing the Petterson norm of Delta:");
Norm = sum(m=0,10,sum(n=0,10-m,Ind(m,n)*(-1)^m*binomial(10,m+n)*binomial(m+n,m)\
                        *delta_r(m)*conj(delta_r(n))))/(3*(-2*I)^11);
print("1) Using Haberland's formula:                  ",Norm);

Norm2 = sum(m=0,10,sum(n=0,10-m,Ind(m,n)*(-1)^m*binomial(10,m+n)*binomial(m+n,m)\
                        *delta_r2(m)*conj(delta_r2(n))))/(3*(-2*I)^11);
print("1') Using Haberland's formula and L-functions: ",Norm2);

Norm3 = 225/(2048*I)*delta_r(1)*delta_r(2);
print("2) Using a rationality theorem:                ",Norm3);