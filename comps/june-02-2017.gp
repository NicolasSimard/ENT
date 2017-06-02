/*Comparing the algebraic delta modular form with the discriminant function of
elliptic curves.*/

L(v) = ellperiods(ellinit(v));
elldisc(v) = ellinit(v).disc;
algdelta(z) = (2*Pi*I)^12*delta(z);

v = [3,1];

print("algdelta(L(",v,")):  ",algdelta(L(v)));
print("disc(E(",v,")):      ",elldisc(v));