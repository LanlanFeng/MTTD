function [X] = Fold_Square( X, dim, i )
N=length(dim);
order=[i:N 1:i-1];
X_temp=reshape(X,dim(order));
X = ipermute(X_temp,order);
end



 