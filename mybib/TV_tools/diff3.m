%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the diff. of one 3-order tensor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function diff_x = diff3(x,  sizeD,w)

tenX   = reshape(x, sizeD);

dfx1     = diff(tenX, 1, 1);
dfy1     = diff(tenX, 1, 2);
dfz1     = diff(tenX, 1, 3);

dfx      = zeros(sizeD);
dfy      = zeros(sizeD);
dfz      = zeros(sizeD);
dfx(1:end-1,:,:) = dfx1;
dfx(end,:,:)     =  tenX(1,:,:) - tenX(end,:,:);
dfy(:,1:end-1,:) = dfy1;
dfy(:,end,:)     = tenX(:,1,:) - tenX(:,end,:);
dfz(:,:,1:end-1) = dfz1;
dfz(:,:,end)     = tenX(:,:,1) - tenX(:,:,end);

diff_x = [w(1)*dfx(:); w(2)*dfy(:); w(3)*dfz(:)];
end
