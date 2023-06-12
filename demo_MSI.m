clc;
clear all;
close all;
addpath(genpath('.\Test_data'));
addpath(genpath('.\mybib'));
filename = ['results_MSI_imshow'];
mkdir(filename);
for i = 2:1:2
    switch i
        case 1
            eval( [ 'load T' ,num2str(13) ]);
            SR = 0.01;
            data_name = 'T13_0.01';
        case 2
            eval( [ 'load T' ,num2str(25) ]);
            SR = 0.03;
            data_name = 'T25_0.03';
    end
    Nway =size(T);
    N = numel(Nway);
    %% Generate known data
    mr = (1-SR)*100;                % Missing ratio (mr);
    Y=zeros(Nway);
    P = round(SR*prod(Nway));
    Omega = randsample(prod(Nway),P);
    Y(Omega)=T(Omega);
    
    j=0;
    %% Method 1： MTTD (FFT)
    j=j+1;
    % choose the parameters
    opts=[];
    alpha=[0.001;1;1]; % for color image
    % alpha=[1;1;1;1]; % for color video
    alphasum = sum( alpha);
    alpha = alpha./ alphasum;
    opts.alpha = alpha;
    f = 0.001;
    opts.gamma =f*alpha;
    opts.rho = 1.1;
    opts.maxIter = 500;
    opts.epsilon = 1e-3;
    opts.T0 = T; % original image
    % Choose transform
    transform.L = @fft; transform.l = 1; transform.inverseL = @ifft;
    % transform.L = @dct; transform.l = 1; transform.inverseL = @idct;
    % main loop
    tic;
    [Xall{j},err] = MTTD_LRTC(Y, Omega, opts, transform);
    % Computing the indexes
    Time(j)=toc;
    [PSNR(j),RSE(j),SSIM(j)]=quality_ll(Xall{j}.*255,T.*255);
    %% Method 2： MTTD (DCT)
    j=j+1;
    % Choose transform
    transform.L = @dct; transform.l = 1; transform.inverseL = @idct;
    tic;
    [Xall{j},err] = MTTD_LRTC(Y, Omega, opts, transform);
    % Computing the indexes
    Time(j)=toc;
    [PSNR(j),RSE(j),SSIM(j)]=quality_ll(Xall{j}.*255,T.*255);
    %% Method 3： MTTD (DCT-Data)
    j=j+1;
    tic;
    [Xall{j},err] = MTTD_LRTCU(Xall{j-1}, Y, Omega, opts, transform);%
    % Computing the indexes
    Time(j)=toc+ Time(j-1);
    [PSNR(j),RSE(j),SSIM(j)]=quality_ll(Xall{j}.*255,T.*255);
    %% Method 4： SMTTD (DCT)
    j=j+1;
    opts=[];
    alpha=[0.001,1,1];
    alphasum = sum( alpha);
    alpha = alpha./ alphasum;
    opts.alpha=alpha;
    opts.gamma =10^-4*alpha;
    opts.beta = 10^-3;
    opts.tau=[1.1,1];
    opts.max_iter=500;
    opts.tol=1e-3;
    opts.T0=T;
    opts.Omega = Omega;
    opts.astv=1;
    opts.lambda=10;
    opts.weight=[1,1,0.5];
    transform.L = @dct; transform.l = 1; transform.inverseL = @idct;
    tic;
    [Xall{j}] = SMTTD_LRTC(Y,opts,transform);
    Time(j)=toc;
    [PSNR(j),RSE(j),SSIM(j)]=quality_ll(Xall{j}.*255,T.*255);
    %% Method 5： MTTD (DCT-Data)
    j=j+1;
    tic;
    [Xall{j},err] = SMTTD_LRTCU(Xall{j-1},Y,opts,transform);
    Time(j)=toc+ Time(j-1);
    [PSNR(j),RSE(j),SSIM(j)]=quality_ll(Xall{j}.*255,T.*255);
  
    save([filename,'\',data_name,'.mat'],'T','Y','Xall','Omega');
end
