quadsetup(-23,2); \\ Define K, pipdata(K), Om_K = CSperiod(K.disc)

/* Petersson norm of \theta_{O_K, ell},
algebrized by the Chowla-Selberg period*/
N1(ell) = pip(pipdata, ell, 1, 1) / Om_K^(4 * ell);

/* Petersson norm of \theta_{O_K, ell},
algebrized by a period Om such that C/Om*O_K is def over H.*/
N2(ell) = pip(pipdata, ell, 1, 1)*(2 * I * Pi * canperiod(K, 1))^(4 * ell);