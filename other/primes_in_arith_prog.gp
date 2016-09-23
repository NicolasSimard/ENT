f(N,M=50000) ={
    my(v=vector(N),piM);
    forprime(p=2,M,if(N%p != 0, v[p%N]+=1));
    piM=sum(n=1,#v,v[n]);
    print("\nThere are ",piM," primes <",M," and "eulerphi(N)," congruence classes.\n");
    printf("Partial distribution: ");
    for(n=1,#v,if(v[n] != 0, printf("%.1f | ", v[n]/piM*100.0)));
    printf("\nUniform distribution: %.1f \n\n",100./eulerphi(N));
}