seq_matrix(out,F,w,p,N,mint,maxt) =
{
    my(L=[],M);
    for(t=mint,maxt,
        print("t=",t);
        L=concat(L,find_seq(F(t),w(t),p,N,250));
    );
    M=matrix(N,maxt-mint+1,n,t,L[N*(t-1)+n]);
    if(out == "pari", return(M));
    out = concat(out,Str("_",p,"_",N,"_",mint,"_",maxt,".md"));
    write("out.tmp",M);
    system(Str("python mat_to_md.py ","out.tmp ",p," ",mint," ",maxt," ",out));
    system("del out.tmp");
}
addhelp(seq_matrix,"seq_matrix(out,F,w,p,N,mint,maxt): Unless out=pari, this "\
"function does not return anything. What it does is compute a matrix M where "\
"M[i,j]=find_min_k(F[t],w[t],p^i) for 1<=i<=N and mint<=t<=maxt. Then "\
"it either prints this matrix to the screen, or calls the python script "\
"mat_to_md.py to convert this matrix to a markdown table and then prints it to out.");
