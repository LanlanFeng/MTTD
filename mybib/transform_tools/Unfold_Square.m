function [X] = Unfold_Square( X, dim, i )
N=ndims(X);
order=[i:N 1:i-1];
X_temp = permute(X,order);
% L=ceil((N-1)/2);
% L=floor((N-1)/2);
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
    elseif dim(3)==2
        %% For 400*200*2  2^4*5^2  2^3*5^2 2
        if i==2
            X=reshape(X_temp,20,20,400);
        elseif i==3
            X=reshape(X_temp,32,25,200);
        else
            X=X_temp;
        end
        
    else
        % %% For 256*256*30
        % % if i==2
        % % X=reshape(X_temp,64,120,256);
        % % elseif i==3
        % %     X=reshape(X_temp,120,64,256);
        % % else
        % %     X=X_temp;
        % % end
        %
        %% For MSI data 512*512*31
        if i==2
            X=reshape(X_temp,128,124,512);
        elseif i==3
            X=reshape(X_temp,124,128,512);
        else
            X=X_temp;
        end
    end
elseif N==4
    if dim(3)==3
        if dim(1)==256
            %%cloud 256*256*3*8
            if i==1
                X=reshape(X_temp,256,256,3*8);
            elseif i==2
                X=reshape(X_temp,32,24,8*256);
            elseif i==3
                X=reshape(X_temp,6,4,256*256);
            else
                X=reshape(X_temp,32,64,256*3);
            end
            
        else
            % suzie_qcif 144*176*3*30
            if i==1
                X=reshape(X_temp,144,176,3*30);
            elseif i==2
                X=reshape(X_temp,22,24,30*144);
            elseif i==3
                X=reshape(X_temp,9,10,144*176);
            else
                X=reshape(X_temp,60,72,176*3);
            end
        end
        
        
        %% For Light Field Data 434*625*9*9
    elseif dim(3)==9
        if i==1
            X = reshape(X_temp,434, 625, 9*9);
        elseif i==2
            X = reshape(X_temp,125, 45, 9*434);
        elseif i==3
            X = reshape(X_temp,9,9,434*625);
        else
            X = reshape(X_temp,63,62,9*625);
        end
    else
        %%cloud 256*256*13*4
        if i==1
            X=reshape(X_temp,256,256,13*4);
        elseif i==2
            X=reshape(X_temp,16*4,13*4,4*256);
        elseif i==3
            X=reshape(X_temp,13,4,256*256);
        else
            X=reshape(X_temp,32,32,256*13);
        end
        
    end
    
end


