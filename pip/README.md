# Petersson inner product of modular forms

## Introduction
This folder contains many formulas to compute the Petersson inner product of the three modular forms delta, delta_11 and delta_5. Those formula come from Cohen's article on Haberland's formula.

For example, the most efficient formula to compute the Petersson norm of delta is the one using its symmetric square L-function. To call it, run GP/PARI in the ENT repo and execute the following command: '\r pip/delta/lsym2.gp'.

## List of formulas implemented
Here is a detailed list of the formulas implemented in each folder

### delta/ folder
- `cohenthm5.1.gp`: Petersson norm of delta using the formula in Thm 5.1 of Cohen's article on Haberland's formula.
- `def.gp`: Petersson norm of delta from the definition as a double integral.
- `haberland.gp`: Petersson norm of delta using Haberlnd's formula.
- `lsym2.gp`: Petersson norm of delta using a special value of the symmetric square L-function of delta (most efficient method so far).
- `periods.gp`: Computation of the periods of delta.
- `rankinl.gp`: Petersson norm of delta as the special value of the Rankin-Selberg convolution L-function of delta with itself.

### delta5/ folder
- `cohenthm5.1.gp`: Petersson norm of delta_5 using the formula in Thm 5.1 of Cohen's article on Haberland's formula.
- `def.gp`: Petersson norm of delta_5 from the definition as a double integral.
- `lsym2.gp`: Petersson norm of delta using a special value of the symmetric square L-function of delta_5 (most efficient method so far).

### delta11/ folder
- `cohenthm5.1.gp`: Petersson norm of delta_11 using the formula in Thm 5.1 of Cohen's article on Haberland's formula.
- `def.gp`: Petersson norm of delta_11 from the definition as a double integral.
- `lsym2.gp`: Petersson norm of delta using a special value of the symmetric square L-function of delta_11 (most efficient method so far).
