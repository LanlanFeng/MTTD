%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the diff. of one 3-order tensor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function diff_x = diff5(x,  sizeD,w)

tenX   = reshape(x, sizeD);

dfx1     = diff(tenX, 1, 1);
dfy1     = diff(tenX, 1, 2);
dfz1     = diff(tenX, 1, 3);
dfl1     = diff(tenX, 1, 4);
dfm1     = diff(tenX, 1, 5);

dfx      = zeros(sizeD);
dfy      = zeros(sizeD);
dfz      = zeros(sizeD);
dfl      = zeros(sizeD);
dfm      = zeros(sizeD);
dfx(1:end-1,:,:,:,:) = dfx1;
dfx(end,:,:,:,:)     =  tenX(1,:,:,:,:) - tenX(end,:,:,:,:);
dfy(:,1:end-1,:,:,:) = dfy1;
dfy(:,end,:,:,:)     = tenX(:,1,:,:,:) - tenX(:,end,:,:,:);
dfz(:,:,1:end-1,:,:) = dfz1;
dfz(:,:,end,:,:)     = tenX(:,:,1,:,:) - tenX(:,:,end,:,:);
dfl(:,:,:,1:end-1,:) = dfl1;
dfl(:,:,:,end,:)     = tenX(:,:,:,1,:) - tenX(:,:,:,end,:);
dfm(:,:,:,:,1:end-1) = dfm1;
dfm(:,:,:,:,end)     = tenX(:,:,:,:,1) - tenX(:,:,:,:,end);
% diff_x=dfx+dfy+dft;
diff_x = [w(1)*dfx(:); w(2)*dfy(:);w(3)*dfz(:);w(4)*dfl(:);w(5)*dfm(:)];
end
