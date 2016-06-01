{
    addhelp(Numth,"This package defines a few number theoretic functions.
    
    *Arithmetic funcitons:
    - ramanujantau(n)->tau(n)
    ");
}

ramanujantau(n) = (5*sigma(n,3)+7*sigma(n,5))*n/12-35*sum(k=1,n-1,(6*k-4*(n-k))*sigma(k,3)*sigma(n-k,5));
