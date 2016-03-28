dtfp(p,N,f) =
{
    my(M=[]);
    for(t=0,10,
        my(dfp);
        dfp=d(bracket(p)(f),-1-t);
        M=concat(M,find_seq(dfp,10-2*t,p,N,250));
    );
    matrix(N,11,n,t,M[N*(t-1)+n])
}

seq_matrix(out,p,i,j,F,w) =
{
    my(M=[],cmd);
    for(t=1,j,
        M=concat(M,find_seq(F(i,j),w(i,j),p,i,250));
    );
    M=matrix(i,j,x,y,M[i*y+x]);
    write("tmp",M);
    cmd = concat("python mat_to_md.py ",p);
    cmd = concat(cmd," tmp");
    system(cmd);
}
