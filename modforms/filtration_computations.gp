seq_matrix(out,F,w,p,N,range) =
{
    my(L=[],M);
    for(t=range[1],range[2],
        L=concat(L,find_seq(F(t),w(t),p,N,250));
    );
    M=matrix(N,range[2]-range[1]+1,n,t,L[N*(t-1)+n]);
    if(out == "pari", return(M));
    write("out.tmp",M);
    system(Str("python mat_to_md.py ","out.tmp ",p," ",range[1]," ",range[2]," ",out));
    system("del out.tmp")
}
addhelp(seq_matrix,"seq_matrix(out,F,w,p,N,range): Unless out=pari, this "\
"function does not return anything. What it does is compute a matrix M where "\
"M[i,j]=find_min_k(F[j],w[j],p^i) for 1<=i<=N and range[1]<=j<=range[2]. Then "\
"it either prints this matrix to the screen, or calls the python script "\
"mat_to_md.py to convert this matrix to a markdown table and then prints it to out.");
