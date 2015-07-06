"""Print a list of candidate discriminants for computations with theta fcts.

The goal is to investigate the Petersson inner product of theta functions
associated with ideal classes in imaginary quadratic fields of prime
discriminant and of small class number. Given a bound N, this small script
prints a list of tuples of the form (-p, h(-p), [f1,...,f_(h(p-))]), where
p < N.
"""

import primes, quadratic

def candidates_list(N):
	L = []
	for p in [p for p in primes.prime_list(N+1) if p % 4 == 3]:
		L.append((-p,quadratic.class_number(-p),
			     quadratic.QuadraticForm.Gaussian_classes_reps(-p)))
	return L

if __name__ == "__main__":
	import sys
	if len(sys.argv) == 1:
		N = int(input("Enter a bound: "))
	else:
		N = int(sys.argv[1])
	for t in candidates_list(N):
		print(t)