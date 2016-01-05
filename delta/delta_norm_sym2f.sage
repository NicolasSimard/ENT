R=RealField(100)
y=R(euler_gamma) #gives the constant with precission 100
print "computing C"
C=2*math.pow(pi,3/2)
k=12

print "computing gammak"
gammak=math.pow(C,-k)*factorial(k-1)
gammakm1=math.pow(C,-k+1)*factorial(k-2)*sqrt(pi)

print "computing bounds"
HnBound=500

HnList=[0 for i in range(HnBound)]
#TESTED
def H(n):
    global HnList
    if n<=0: return 0
    if HnList[n]!=0:
        return HnList[n]
    HnList[n]=HnList[n-1]+1/n
    return HnList[n]

# FactList=[1 for i in range(FactBound)]
# THERE SEEMS TO BE A PROBLEM...
# def fact(n):
#    if n==0: 
#        return 1
#    global FactList
#    if FactList[n]!=1:
#        return FactList[n]
#    FactList[n]=FactList[n-1]*n
#    return FactList[n]

print "computing delta"
delta=CuspForms(Gamma0(1),12).basis()[0]
print "computing coefficients"
a=delta.coefficients(range(100))

def Omega(n):
    return sum(d[1] for d in list(factor(n)))

def A(n):
    return sum(a[int(n/m)]**2*m**(k-1)*(-1)**Omega(m) for m in divisors(n))

# We add a 0 at the begining of the list so that the indices match
# Note that the bound is (k-2)/2+1=k/2 because range drops the last integer
# F1kn=[0]
# F1kn=F1kn+[(-1)**(m+1)*factorial(2*m-1)*math.pow(C,-2*m)/factorial(int(k/2)-m-1) for m in range(1,int(k/2))]
# WORKS (BUT NOT TESTED)
def F1k(s,x):
    return sum((-1)**(m+1)*factorial(2*m-1)*math.pow(C*x,-2*m)/(factorial(int(k/2)-m-1)*(s-2*m)) for m in range(1,int(k/2)))

def F2k_coeff(s,x,m):
    return (-1)**(m+1)*2**(2*m+k)*factorial(m+int(k/2))*math.pow(C*x,2*m+1)/(factorial(2*m+1)*factorial(2*m+k)*(s+2*m+1))

# F2kn=[(-1)**(m+1)*2**(2*m+k)*factorial(m+int(k/2))*math.pow(C,2*m+1)/(factorial(2*m+1)*factorial(2*m+k)) for m in range(0,M2)]
# WORKS (BUT NOT TESTED)
def F2k(s,x):
    S=0
    m=0
    while abs(F2k_coeff(s,x,m))>Err:
        S+=F2k_coeff(s,x,m)
        m+=1
    return S

def F3k_coeff(s,x,m):
    return (-1)**(m+1)*math.pow(C*x,2*m)/(factorial(2*m)*factorial(m+int(k/2)-1)*(s+2*m))*(2*H(2*m)+H(m+int(k/2)-1)-3*y-2*log(C*x)+2/(2*m+s)) 

# F3kn=[(-1)**(m+1)*math.pow(C,2*m)/(factorial(2*m)*factorial(m+int(k/2)-1)) for m in range(0,M3)]
def F3k(s,x):
    # L=[2*H(2*m)+H(m+int(k/2)-1)-3*y-2*log(C*x)+2/(2*m+s) for m in range(0,M3)]
    S=0
    m=0
    while abs(F3k_coeff(s,x,m))>Err:
        S+=F3k_coeff(s,x,m)
        m+=1
    return S

def Fkk(x): # Formula of Fk(s,x) with s=k
    return gammak-x**k*(2*F1k(k,x)+sqrt(pi)*F2k(k,x)+F3k(k,x))

def Fkkm1(x): # Formula of Fk(s,x) with s=k-1
    return gammakm1-x**(k-1)*(2*F1k(k-1,x)+sqrt(pi)*F2k(k-1,x)+F3k(k-1,x))

print"C=",C,"; k=",k,"; gamma(k)=",gammak
# for x in range(M):
#    m=0
#    while abs(F2k_coeff(k,x,m))>Err: m+=1
#    print "x=",x,": m=",m
# print m,": ",F2k_coeff(k,2,m)
print"Computing the petersson norm... "
Err=1e-2
S=0
n=1
while N(abs(A(n)/(n**k)*(Fkk(n)+n*Fkkm1(n))))>Err:
    print "term ",n,": ",N(A(n)/(n**k)*(Fkk(n)+n*Fkkm1(n)),20)
    S+=A(n)/(n**k)*(Fkk(n)+n*Fkkm1(n))
    n+=1
    print 2**(1-k)*math.pow(pi,5)*S
Peter=2**(1-k)*math.pow(pi,5)*S
print Peter
