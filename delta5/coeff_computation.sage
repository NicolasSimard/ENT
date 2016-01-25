import json

def prime_squares(N):
    return [p^2 for p in prime_range(int(sqrt(N)))]

delta5 = Newforms(Gamma0(5),4)[0]

def compute_coeffs(bound):
    L = delta5.coefficients(prime_squares(bound))
    coeff = [0]*bound
    count = 0
    for p in prime_range(int(sqrt(bound))):
        coeff[p^2-1] = L[count]
        count += 1

    coeff_str = "a = " + str(coeff) + ";"
    json.dump(coeff_str,open('partial_coeff.data','w'))
