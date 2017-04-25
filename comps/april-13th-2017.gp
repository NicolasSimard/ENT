/*Experiment with Lang's statements in Chapter 12 of his book on elliptic functions.*/

phi_alpha(alpha,z) = matdet(alpha)^12*delta([alpha[1,1]*z+alpha[1,2],alpha[2,1]*z+alpha[2,2]])/delta(z);

alpha = [3,5;1,3];
z = sqrt(-5);

x=phi_alpha(alpha,z);
print("Computing algdep with degree 30...");
algdep(x,30)