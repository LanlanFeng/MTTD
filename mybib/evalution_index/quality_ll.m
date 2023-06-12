function [PSNR, RSE, SSIM] = quality_ll(imagery1, imagery2)
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
PSNR = zeros(Nway(3),1);
RSE = PSNR;
SSIM = PSNR;
for i = 1:Nway(3)
    PSNR(i) = psnr(uint8(imagery1(:, :, i)), uint8(imagery2(:, :, i)));
%     PSNR(i) = psnr(imagery1(:, :, i), imagery2(:, :, i));
%     RSE(i) = norm(imagery1(:, :, i)-imagery2(:, :, i),'fro')/norm(imagery2(:, :, i),'fro');
    SSIM(i) = ssim(uint8(imagery1(:, :, i)), uint8(imagery2(:, :, i)));%dynamic range is 0-255
end
PSNR = mean(PSNR);
% RSE = mean(RSE);
RSE = norm(imagery1(:)-imagery2(:),'fro')/norm(imagery2(:),'fro');
SSIM = mean(SSIM);

