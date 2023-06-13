function [X] = Unfold_MTTD_RTPCA( X, dim, i )
N=ndims(X);
order=[i:N 1:i-1];
% L=ceil((N-1)/2);
% L=floor((N-1)/2);
X_temp = permute(X,order);
% X=reshape(X_temp,prod(dim(order(1:L))),prod(dim(order(L+1:N-1))),dim(order(N)));%TR
X=reshape(X_temp,dim(order(1)),dim(order(2)),prod(dim(order(3:N))));%TUCKER
%% For color image
if N==3
    if dim(3)==3
        %% For 256*256*3
        if i==2
            X=reshape(X_temp,32,24,256);
        elseif i==3
            X=reshape(X_temp,24,32,256);
        else
            X=X_temp;
        end
        
    else
        
        %% For 120*160*30   2^3*3*5 2^5*5 3 2*3*5
        if i==2
            X=reshape(X_temp,2^4*5,2^2*3*5,120);
        elseif i==3
            X=reshape(X_temp,2^2*3*5,2^2*3*5, 160);
        else
            X=X_temp;
        end
        
    end
end
end


