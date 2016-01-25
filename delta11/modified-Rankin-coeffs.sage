attach('coeff_computation.sage')
attach('convolution.sage')

def ind(n):
    if gcd(11,n) == 1 and int(sqrt(n))^2 == n:
        return 1
    return 0

def compute_Rankin_coeffs(bound,file_name = ""):
    convol_coeffs = [a*conjugate(a)/(n+1) for (n,a) in enumerate(compute_coeffs(bound))]
    coeff = dir_convolution([0]+[ind(n) for n in range(1,bound)],[0]+convol_coeffs)[1:]
    if len(file_name) == 0:
        return coeff
    else:
        coeff_str = "a = " + str(coeff) + ";"
        json.dump(coeff_str,open(file_name,'w'))

