clc;
clear all;
close all;
addpath(genpath('.\Test_data'));
addpath(genpath('.\mybib'));
filename = ['results_colorimage_denoising'];
mkdir(filename);

%% Load initial data
load('lena.mat');
data_name = 'lena';
X = T;
if max(X(:))>1
    X = my_normalized(X);
end

Nway = size(X);
Ndim = ndims(X);

rhos = 0.7;  %% The corrupted ratio
fprintf('=== The sample ratio is %4.2f ===\n', rhos);
Xn = imnoise(X,'salt & pepper',rhos);


%% Method: MTTD
i=0;
i = i+1;

% set the parameters
opts=[];
if rhos==0.1
    alpha=[1, 1, 1];
    opts.alpha=alpha/sum(alpha(:));
    opts.tol = 3e-3 ; 
end
if rhos==0.3
    alpha=[1, 1, 1];
    opts.alpha=alpha/sum(alpha(:));
    opts.tol = 5e-4 ;
end
if rhos==0.5
    alpha=[1, 0.5, 0.5];
    opts.alpha=alpha/sum(alpha(:));
    opts.tol = 2e-3;
end
if rhos==0.7
    alpha=[1, 0.5, 0.5];
    opts.alpha=alpha/sum(alpha(:));
    opts.tol = 1e-4 ;
end
opts.maxIter=5000;
opts.rho = [1,1.0];
omega = 500;
opts.gamma = opts.alpha/omega;
opts.beta = 1/100;
opts.max_gamma = 1e10*ones(Ndim);
opts.lambda = Set_lambdall(X,Nway,opts.alpha);
opts.T0=X;
%% Choose transform
transform.L = @fft; transform.l = 1; transform.inverseL = @ifft;
% transform.L = @dct; transform.l = 1; transform.inverseL = @idct;
t0= tic;
[Xall{i},~,~]=MTTD_RTPCA(Xn,opts,transform);
L = abs(Xall{i});
Xall{i} = L;
figure,imshow(Xall{i});
Time(i)= toc(t0);
[PSNRA(i), RSEA(i), SSIMA(i)] = quality_ll(Xall{i}*255, X*255);



savePath= [filename,'\', data_name,'rhos', num2str(rhos),'.mat'];
save(savePath, 'PSNRA','RSEA','SSIMA','Time','Xall','Xn','X');  
