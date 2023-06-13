function lambda = Set_lambdall(T,Nway,alpha) 
dim = size(T);
lambda = 0;
for i=1:ndims(T)
    OutX = Unfold_MTTD_RTPCA( T, dim, i );
    [N1,N2,N3]=size(OutX);
    temp=alpha(i)/sqrt(max(N1,N2)*N3);
    lambda = lambda + temp;
end