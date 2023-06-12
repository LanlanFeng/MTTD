%% test 8 color images 
clc;
clear all;
close all;
addpath(genpath('.\Test_data'));
addpath(genpath('.\mybib'));
filemame = ['results_colorimage'];
mkdir(filemame);
for data_num =1:1:8
    switch data_num
        case 1
            load('airplane.mat');
            data_name = 'airplane';
        case 2
            load('sails.mat');
            data_name = 'sails';
        case 3
            load('barbara.mat');
            data_name = 'barbara';
        case 4
            load('house.mat');
            data_name = 'house';
        case 5
            load('lena.mat');
            data_name = 'lena';
        case 6
            load('peppers.mat');
            data_name = 'peppers';
        case 7
            load('sailboat.mat');
            data_name = 'sailboat';
        case 8
            load('tulips.mat');
            data_name = 'tulips';
    end
    Nway =size(T);
    N = numel(Nway);
    for i=1:1:1
        SR =i/20;                    % Sample ratio (SR), e.g. 0.1 = 10% known samples
        mr = (1-SR)*100;             % Missing ratio (mr);
        Y=zeros(Nway);
        % Generate obsvered data
        P = round(SR*prod(Nway));
        Omega = randsample(prod(Nway),P);
        Y(Omega)=T(Omega);
        %         O = zeros(Nway);
        %         O(Omega) = 1;
        j=0;
        %% Method 1： MTTD (FFT)
        j=j+1;
        % choose the parameters
        opts=[];
        alpha=[0.01;1;1]; % for color image
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
        Time(j,i)=toc;
        [PSNR(j,i),RSE(j,i),SSIM(j,i)]=quality_ll(Xall{j}.*255,T.*255);
        %% Method 2： MTTD (DCT)
        j=j+1;
        % Choose transform
        transform.L = @dct; transform.l = 1; transform.inverseL = @idct;
        tic;
        [Xall{j},err] = MTTD_LRTC(Y, Omega, opts, transform);
        % Computing the indexes
        Time(j,i)=toc;
        [PSNR(j,i),RSE(j,i),SSIM(j,i)]=quality_ll(Xall{j}.*255,T.*255);
        %% Method 3： MTTD (DCT-Data)
        j=j+1;
        tic;
        [Xall{j},err] = MTTD_LRTCU(Xall{j-1}, Y, Omega, opts, transform);%
        % Computing the indexes
        Time(j,i)=toc+ Time(j-1,i);
        [PSNR(j,i),RSE(j,i),SSIM(j,i)]=quality_ll(Xall{j}.*255,T.*255);
        %% Method 4： SMTTD (DCT)
        j=j+1;
        opts=[];
        alpha=[0.01,1,1];
        alphasum = sum( alpha);
        alpha = alpha./ alphasum;
        opts.alpha=alpha;
        opts.gamma =10^-4*alpha;
        opts.beta = 10^-2;
        opts.tau=[1.1,1];
        opts.max_iter=500;
        opts.tol=1e-3;
        opts.T0=T;
        opts.Omega = Omega;
        opts.astv=1;
        opts.lambda=1;
        opts.weight=[1,1,0];
        transform.L = @dct; transform.l = 1; transform.inverseL = @idct;
        tic;
        [Xall{j}] = SMTTD_LRTC(Y,opts,transform);
        Time(j,i)=toc;
        [PSNR(j,i),RSE(j,i),SSIM(j,i)]=quality_ll(Xall{j}.*255,T.*255);
        %% Method 5： MTTD (DCT-Data)
        j=j+1;
%         opts.weight=[0.5,0.5,0];
%         opts.weight=[1,1,0];
%         opts.tol=2e-3;        
        tic;
        [Xall{j},err] = SMTTD_LRTCU(Xall{j-1},Y,opts,transform);
        Time(j,i)=toc+ Time(j-1,i);
        [PSNR(j,i),RSE(j,i),SSIM(j,i)]=quality_ll(Xall{j}.*255,T.*255);
    end
    savePath=[filemame,'\', data_name,'_varySR','j',num2str(j),'.mat'];
    save(savePath, 'PSNR','RSE','SSIM','Time');  %save as mat file
end
