/* A computation session with the theta series of weight 1 attached to the
quadratic field of discriminant -23. To see the computations, open PARI/GP and
type ?> \r session12_02_16.gp. Do not use ?> read("session12_02_16.gp")
because you won't see the output!

Scripts needed:
1) ENT/theta/theta_function.gp
2) ENT/quadratic.gp
*/

read("theta_function.gp");

/* This function was taken in the computel script.*/
errprint(x)={
    if(type(x)=="t_COMPLEX",x=abs(x));
    if(x==0,return(concat("1E-",default(realprecision)+1)),
    return(concat(concat(truncate(x/10^floor(log(abs(x))/log(10))),"E"),floor(log(abs(x))/log(10)))));
}

\p 200
allocatemem(32000000);
print("*** Setup ***");
print("-- Let D=",D=-23,".");
print("-- The reduced forms are :",forms=reduced_forms(D));
print("-- Define the theta functions F0 and F1 of weight 1 attached to\n",\
      "   ",f0=forms[1]," and "f1=forms[2]);
print(">> F0 = ",theta_qexp(f0,0,15,"series"));
print(">> F1 = ",theta_qexp(f1,0,15,"series"));
F0(z) = theta_function(f0,0,z,100000);
F1(z) = theta_function(f1,0,z,100000);
f2=forms[3];
F2(z) = theta_function(f2,0,z,100000);
print("-- These theta functions have weight 1, level 23 and their\n",\
      "   character is the kronecker symbol (-23/-).");
print("-- The Heegner forms corresponding to the quadratic forms above\n",\
      "   ",HFs=vector(3,n,Heegner_form(forms[n],23,23))," respectively");
print("-- Let t0 and t1 ne the roots of ",HFs[1]," and ",HFs[2]," respectively.")
t0 = tau(HFs[1]);
t1 = tau(HFs[2]);
t2 = tau(HFs[3]);
print("-- Note that the root of ",HFs[2]," is -conj(t1).");

print("\n\n*** Computations ***");
print("-- We now evaluate a0001=F0(t0)/F0(t1). The theory says that it is algebraic.");
a0001 = F0(t0)/F0(t1);
print("-- We now try to find its minimal polynomial and find the irreducible factors.");
p0001 = factor(algdep(a0001,15));
print(">> ",p0001[,1]~);
print("-- Then we evaluate each factor at a0001. This gives:");
print(">> ",apply(p->errprint(subst(p,x,a0001)),p0001[,1]~));
print("-- So the minimal polynomial should be: \n",\
      ">> ",m0001=p0001[2,1]);

print("\n")
print("-- We now evaluate a0010=F0(t0)/F1(t0). Again, we know it is algebraic.");
a0010 = F0(t0)/F1(t0);
print("-- We now try to find the minimal polynomial of a0010 and we factor it.");
p0010 = factor(algdep(a0010,15));
print(">> ",p0010[,1]~);
print("-- Then we evaluate each factor at a0010. This gives:");
print(">> ",apply(p->errprint(subst(p,x,a0010)),p0010[,1]~));
print("-- So the minimal polynomial should be: \n",\
      ">> ",m0010=p0010[2,1]);

print("\n")
print("-- We now evaluate a0011=F0(t0)/F1(t1). Again, we know it is algebraic.");
a0011 = F0(t0)/F1(t1);
print("-- We now try to find its minimal polynomial and print the irreducible factors.");
p0011 = factor(algdep(a0011,15));
print(">> ",p0011[,1]~);
print("-- Then we evaluate each factor at a0011. This gives:");
print(">> ",apply(p->errprint(subst(p,x,a0011)),p0011[,1]~));
print("-- So the minimal polynomial should be: \n",\
      ">> ",m0011=p0011[4,1]);

print("\n")
print("-- We now evaluate a1011=F1(t0)/F1(t1). Again, we know it is algebraic.");
a1011 = F1(t0)/F1(t1);
print("-- We now try to find its minimal polynomial and print the irreducible factors.");
p1011 = factor(algdep(a1011,10));
print(">> ",p1011[,1]~);
print("-- Then we evaluate each factor at a1011. This gives:");
print(">> ",apply(p->errprint(subst(p,x,a1011)),p1011[,1]~));
print("-- So the minimal polynomial should be: \n",\
      ">> ",m1011=p1011[3,1]);

