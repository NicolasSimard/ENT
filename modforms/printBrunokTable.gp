table(p,N,j) =
{
    printf("p^n\\j  ");
    for(i=0,j,printf("|%-7i",i));
    printf("\n");
    printf("-------");
    for(i=0,j,printf("|-------"));
    printf("\n");
    for(n=1,N,
        printf("%i^%i    ",p,n);
        for(i=0,j,
            printf("|%-7i",find_min_k(V(p,i)(E2star(p)),2,p^n)[1]); 
        );
        printf("\n");
    );
    
}