# Theta series and their Petersson inner product.

## Introduction
This is where most computations with theta functions take place. The main
script is `ENT/theta/Thetapip.gp`, which defines a bunch of functions to efficiently compute
the Petersson inner product of binary theta series attached to imaginary quadratic
fields. A detailed presentation of the formulas that are used in Thetapip.gp can
be found in the `Notes/` repo.

## Content of `ENT/theta/` directory
Here are some of the most useful scripts in the directory.

- `ENT/theta/Thetapip.gp`: Defines a bunch of functions to compute the Petersson inner
product of theta series attached to ideals of Hecke characters of imaginary
quadratic fields.
- `ENT/theta/Thetafunc.gp`: Defines the theta series themselves (their q-expansion, etc.)

## Examples
Here are some examples of computations one can do with the scripts.

### Computing Petersson inner products of theta series attached to ideals
Here is an example of a computations with theta series attached to
ideals in the field $K=\Q(\sqrt(-23))$.

#### Initialize
First of all, start PARI/GP in the `ENT/theta/` directory. Then load the script,
set the working precision and define K.

```
(14:16) gp > \r Thetapip.gp
(14:16) gp > \p 500
   realprecision = 500 significant digits
(14:16) gp > K=bnfinit('x^2+23);
```
Before we start computing, we need to initialize the Petersson inner product by
calling `pipinit(K)`.

```
(14:16) gp > data=pipinit(K);
```

#### Start computations
Then we can, for example, compute the Petersson norm of the theta series
attached to the ideal 1 for ell=2:

```
(14:17) gp > ell=2;
(14:17) gp > pip(data ,ell,1,1)
%171 = 0.0043280662913848413894386[...] + 0.E-505*I
```

where the `[...]` means that I cut the output to make it fit on this page (the
PARI output is much longer)! Of course, this becomes more interesting if we
normalize by the Chowla-Selberg period (defined in the `ENT/Quadratic.gp` script),
since we know that we obtain algebraic integers:

```
(14:17) gp > \r ../Quadratic.gp
(14:17) gp > OmK=CSperiod(-23);
(14:17) gp > algdep(pip(data,ell,1,1)/OmK^(4*ell),10)
%205 = x^9 - 61143680*x^6 + 793604147138560*x^3 - 3118852808507803303936
```

We can then see that the cube of the roots of this polynomial generate the
Hilbert class field of K:

```
(14:17) gp > polredbest(%205)
%213 = x^9 - x^3 - 1
```

Here is an example with theta series attached to non-trivial ideals:

```
(14:30) gp > p2=idealprimedec(K,2)[1];
(14:31) gp > p3=idealprimedec(K,3)[1];
(14:31) gp > pip(data,ell,p2,p3)
%227 = 0.0879462923155238887304235[...] - 0.0116257152983884058845420[...]*I
```

That's it!

### Computing Petersson inner products of theta series attached to Hecke characters
Here is an example of a computations with theta series attached to Hecke characters
of the field $K=\Q(\sqrt(-23))$.

#### Initialize
As in the previous example, start PARI/GP in the `ENT/theta/` directory, load the scripts,
set the working precision, define K and call `pipinit(K)`:

```
(14:39) gp > \r Thetapip.gp
(14:39) gp > \r ../Quadratic.gp
(14:39) gp > \r ../lfunc/nf/quadhecke.gp
(14:40) gp > \p 500
   realprecision = 500 significant digits
(14:40) gp > K=bnfinit('x^2+23);
(14:40) gp > data=pipinit(K);
```

Note that we have to load an extra script, namely `ENT/lfunc/nf/quadhecke.gp`, which defines the L-function attached to Hecke characters of imaginary quadratic fields. The script `ENT/Quadratic.gp` is also needed to define Hecke characters.

We can now compute the Petersson norm of the theta series attached to one of
these Hecke characters of infinity type 2*ell for ell=2:

```
(14:40) gp > ell=2;
(14:40) gp > pnorm(data,[[0],[2*ell,0]])
%253 = 0.0016833357054580643163186[...] - 1.845288302 E-503*I
```

Again, we can normalize by the Chowla-Selberg period to obtain algebraic numbers.

```
(14:41) gp > OmK=CSperiod(-23);
(14:41) gp > algdep(pnorm(data,[[0],[2*ell,0]])/OmK^(4*ell),20)
%287 = x^9 - 188333424*x^6 + 639297099696618*x^3 - 389589775351907071941
(14:41) gp > algdep(pnorm(data,[[1],[2*ell,0]])/OmK^(4*ell),20)
%288 = x^9 - 188333424*x^6 + 639297099696618*x^3 - 389589775351907071941
(14:41) gp > algdep(pnorm(data,[[2],[2*ell,0]])/OmK^(4*ell),20)
%289 = x^9 - 188333424*x^6 + 639297099696618*x^3 - 389589775351907071941
```
