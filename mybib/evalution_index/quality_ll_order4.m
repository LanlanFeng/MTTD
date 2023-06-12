function [PSNR, RSE, SSIM] = quality_ll_order4(imagery1, imagery2)
%==========================================================================
% Evaluates the quality assessment indices for two tensors.
%
% Syntax:
%   [psnr, rse, ssim] = quality(imagery1, imagery2)
%
% Input:
%   imagery1 - the target tensor
%   imagery2 - the reference tensor

% NOTE: the tensor is a M*N*K array and DYNAMIC RANGE [0, 255].

% Output:
%   psnr - Peak Signal-to-Noise Ratio
%   rse  - Relative Square Error
%   ssim - Structure Similarity

Nway = size(imagery1);
PSNR = zeros(Nway(3),Nway(4));
RSE = PSNR;
SSIM = PSNR;
for i = 1:Nway(3)
    for j = 1:Nway(4)
        PSNR(i, j) = psnr(uint8(squeeze(imagery1(:, :, i, j))), uint8(squeeze(imagery2(:, :, i, j))));
        %     PSNR(i) = psnr(imagery1(:, :, i), imagery2(:, :, i));
        %     RSE(i) = norm(imagery1(:, :, i)-imagery2(:, :, i),'fro')/norm(imagery2(:, :, i),'fro');
        SSIM(i, j) = ssim(uint8(imagery1(:, :, i, j)), uint8(imagery2(:, :, i, j)));%dynamic range is 0-255
        
    end
end
PSNR = mean(PSNR(:));
% RSE = mean(RSE);
RSE = norm(imagery1(:)-imagery2(:),'fro')/norm(imagery2(:),'fro');
SSIM = mean(SSIM(:));

