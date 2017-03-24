# ENT
A repository which contains PARI/GP code (mainly) to experiment with modular forms
and other number theoretic objects. The code was tested with PARI 2.7.2.

## Introduction
This repo mainly contains scripts to do computations with modular forms. A big
part of the code is focused on computing the Petersson inner product of binary
theta series attached to imaginary quadratic fields (see theta/ directory).
Another project was about p-adic modular forms (with Bruno Joyal, see padic/
directory).

The code is organized in the following way: scripts that start with a capital
letter contain a collection of related functions. For example, the Modforms.gp
script contains functions related to modular forms (Eisenstein series, delta
function, Victor Miller basis, etc.). In general, each function has a help text,
which you can view by typing ?function in PARI/GP. Moreover, each capitalized
script (like Modform.gp, for example) has a general help function, which you can
view by typing `? + 'Name of the script'` (like `?Modform`, for example). This
help text gives a general description of the script, lists all the available
functions and gives their signature. This is very useful if you don't remember
which function returns the Victor Miller Basis in Modform.gp, for example! See
the examples below.

Note that this repo also contains Python 3 code, which I often used to
manipulate text data (to produce Markdown files, for example).

## Content of ENT repository
Here is a short list of some of the available scripts in the repo. Recall that
you can obtain a list of available functions in each capitalized script by
typing `? + 'Name of the script'` in PARI/GP.

### General scripts
- Modform.gp: Functions related to modular forms.
- Quadratic.gp: Functions related to positive definite quadratic forms or
imaginary quadratic fields.
- Jacform.gp: Functions related to Jacobi forms.

### Some projects
Each directory contains a README file with a description of the project.

- theta/ directory: Theta series and their Petersson inner product.
- lfunc/ directory: L-functions of modular forms.
- padic/ directory: filtration of p-adic modular forms.
- pip/ directory: Petersson inner product of some specific modular forms.

## Examples of use
We present here some examples of use of the general scripts of the repo. For
in the different projects, see the README files of those projects.

We illustrate the use of the Modform.gp script.

### Reading the script
Start PARI/GP in the ENT repo and simply use `\r` or `read` commands, like that:
`\r Modform.gp;` or `read('Modform.gp');`.

### General help for the script
After loading the script, type `?Modform` (without the '.gp' extension!) to
obtain general help for the script. For example:

```
(11:41) gp > \r Modform.gp
(11:41) gp > ?Modform
This package defines some common modular forms, some operators
on them and tools to compute basis of certain spaces of modular forms.

In general, a modular form is represented by a q-expansion. The precision
of the q-expansion is determined by the constant default(seriesprecision)
which can be changed with the command ps prec.

*General modular forms:
- G(k,'q) -> q-expansion; G(k,z) -> G_k(z) (numerical value)
- E(k,'q) -> q-expansion; E(k,z) -> E_k(z) (numerical value)
- G2star(z) -> numerical value
- theta0(x) -> q-expansion
- theta1(x) -> q-expansion
- ellj('q) -> q-expansion (PARI built-in); ellj(z) -> j(z) (numerical value)
- delta('q) -> q-expansion; delta(z) -> delta(z) (numerical value)
- fd(d) -> q-expansion
- gD(D) -> q-expansion

*Basis of modular forms:
- vmbasis(k,N,flag) -> [f1('q),f2('q),...,fd('q)] (q-expansions of Victor
Miller basis)
- GktoG4G6(k) -> P(G4,G6) : P(G4,G6)=Gk

*Operators:
- U(n,f('q)) -> U_n(f('q)) (U_n operator in level 1)
- V(n,f('q)) -> V_n(f{'q)) (V_n operator in level 1)
- pstab(p,f('q)) -> f('q)^[p] (p-stabilisation of f)
- rcbracket(f('q),k_f,g('q),k_g) -> [f('q),g('q)] (Rankin-Cohen bracket)
- dop(f('q),n) -> d^n(f('q)) (d='q*d/d'q)
- dopformal(P,n) -> d^n(P) (P in C['E2,'E4,'E6])
- delkformal(P,n) -> del_k^n(P) (P in C['E2s,'E4,'E6, 'G2s,'G4,'G6])

*Other functions:
- jpol(p('q)) -> P(X): P(j('q)) = p('q)
```

### Help for a specific function
Once you know which functions are available, you may want to obtain specific
about a particular function. In this case, simply type `? + 'Name of the function'`:

```
(11:41) gp > ?jpol
jpol(f): Given a q-expansion f(q) with enough coefficients,
returns a polynomial g(Y) such that g(j)=f, where j is the j-function.
```
