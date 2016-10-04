# L-function experiments

## Introduction
This directory contains computations related to L-functions in general. So far, it contains computations related to L-functions of modular forms (in the `ENT/lfunc/mf/` directory) and number fields (in the `ENT/lfunc/nf/` directory).

## Content of the `ENT/lfunc/` directory

- `ENT/lfunc/mf/` directory: L-function of the delta function, symmetric square L-fuction of delta and Rankin-Selberg convolution L-function of a certain class of newforms. Those scripts use the `ENT/computel` script.

- `ENT/lfunc/nf/` directory: Dirichlet L-function, Dedekind zeta function and Hecke L-function of imaginary quadratic fields. Also contains examples. Some of those scripts use the `ENT/computel` script, but some don't.

## A remark on the Hecke L-function of imaginary quadratic fields
The `ENT/lfunc/nf/quadhecke.gp` script is an implementation of the formulas in my notes `Notes/Theta Norm`, while `ENT/lfunc/nf/quadhecke-computel.gp` uses the `ENT/computel` script. For large integral arguments, the formulas in the first script are much more efficient.
