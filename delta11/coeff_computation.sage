import json

def prime_squares(N):
    return [p^2 for p in prime_range(int(sqrt(N)))]

delta11 = Newforms(Gamma0(11),2)[0]

def compute_partial_coeffs(bound,file_name = ""):
    L = delta11.coefficients(prime_squares(bound))
    coeff = [0]*bound
    count = 0
    for p in prime_range(int(sqrt(bound))):
        coeff[p^2-1] = L[count]
        count += 1

    if len(file_name) == 0:
        return coeff
    else:
        coeff_str = "a = " + str(coeff) + ";"
        json.dump(coeff_str,open(file_name,'w'))

def compute_coeffs(bound, file_name = ""):
    coeff = delta11.coefficients(bound)
    if len(file_name) == 0:
        return coeff
    else:
        coeff_str = "a=" + str(coeff) + ";"
        json.dump(coeff_str,open(file_name,'w'))
