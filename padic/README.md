# Filtration of p-adic modular forms
Comput ethe filtration of a p-adic modular form.

## Introduction
This project was initiated by a questionning of Bruno Joyal on p-adic modular 
forms. This lead to a script that tries to compute the filtration of p-adic 
modular forms. 

Given a p-adic modular form a la Serre f (i.e. a q-expansion) of weight 
\kappa and valuation 0, and a positive integer N, it is known that f is 
congruent modulo p^N to a classical modular form of some weight k. Moreover, 
k is congruent to \kappa modulo phi(p^N). Since there are no classical 
modular forms of negative weight, there is a smallest weight k with this 
property. This weight is denoted k_N and the sequence (k_N)_{N>=1} is called 
the filtration of f. It is also convenient to define the normalized 
filtration as (h_N)_{N>=1}, where h_N=(k_N-w)/phi(p^N) (which is an integer...)

## Main scripts in the folder
The most important script is `Filtration.gp`, which computes the filtration 
of a given p-adic modular form. The results of some computations are 
presented in `filtration-data.md`. The other scripts are mostly utilitary 
scripts. For example, the script `mat-to-md.py` is a Python 3 script which 
takes a input a PARI/GP matrix and returns it as a markdown table. It was 
very useful to build the `filtration-data.md` file! 