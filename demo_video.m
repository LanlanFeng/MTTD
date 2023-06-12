clc;
clear all;
close all;
addpath(genpath('.\Test_data'));
addpath(genpath('.\mybib'));
filename = ['results_video'];
mkdir(filename);
for datanum =1:1:2
    switch datanum
        case 1
            load('suzie_qcif.mat')%
            T = im2double(suzie_qcif(:,:,:,41:70));
            data_name = 'suzie_qcif';
        case 2
            load('foreman_qcif.mat')%144*176*3*30
            T = im2double(foreman_qcif(:,:,:,1:30));
            data_name = 'foreman_qcif';
            
    end
    Nway =size(T);
    N = numel(Nway);
    %% Generate Omega data
    for i=1:2:3
        SRall = [0.01,0.05,0.1];
        SR = SRall(i);                    % Sample ratio (SR), e.g. 0.1 = 10% known samples
        mr = (1-SR)*100;                % Missing ratio (mr);
        %         load ([filename,'\', data_name,'_SR',num2str(SR),'a.mat']);
        Y=zeros(Nway);
        % Generate known data
        P = round(SR*prod(Nway));
        Omega = randsample(prod(Nway),P);
        Y(Omega)=T(Omega);
        O = zeros(Nway);
        O(Omega) = 1;
        
        j=0;
        %% Method 1: MTTD (FFT)
        j=j+1;
        % choose the parameters
        opts=[];
        alpha=[1;1;1;1]; % for SR=0.05
        alphasum = sum( alpha);
        alpha = alpha./ alphasum;
        opts.alpha = alpha;
        f = 1e-4;
        opts.gamma =f*alpha;
        opts.rho = 1.1;
        opts.maxIter = 500;
        epsilonall = [1e-5,1e-3,1e-3];
        opts.epsilon = epsilonall(i);
        opts.T0 = T; % original image
        % Choose transform
        transform.L = @fft; transform.l = 1; transform.inverseL = @ifft;
        % main loop
        tic;
        [Xall{j}] = MTTD_LRTC(Y, Omega, opts, transform);
        Time(i,j)=toc;
        [PSNR(i,j),RSE(i,j),SSIM(i,j)]=quality_ll_order4(Xall{j}.*255,T.*255);
        %% Method 2: MTTD (DCT)
        j=j+1
        alpha=[1;1;1;1];
        alphasum = sum( alpha);
        alpha = alpha./ alphasum;
        opts.alpha = alpha;
        f = 1e-4;
        opts.gamma =f*alpha;
        epsilonall = [1e-5,1e-3,1e-3];
        opts.epsilon = epsilonall(i);
        transform.L = @dct; transform.l = 1; transform.inverseL = @idct;
        tic;
        [Xall{j},err] = MTTD_LRTC(Y, Omega, opts, transform);
        Time(i,j)=toc;
        [PSNR(i,j),RSE(i,j),SSIM(i,j)]=quality_ll_order4(Xall{j}.*255,T.*255);
        %% Method 3: MTTD (Data)
        j=j+1
        tic
        if PSNR(i,j-1)>=PSNR(i,j-2)
            U=Xall{j-1};
        else
            U=Xall{j-2};
        end
        epsilonall = [1e-3,2e-3,2e-3];
        [Xall{j},err] = MTTD_LRTCU(U , Y, Omega, opts, transform);
        % 指标计算
        Time(i,j)=toc + Time(i,j-1);
        [PSNR(i,j),RSE(i,j),SSIM(i,j)]=quality_ll_order4(Xall{j}.*255,T.*255);
        %% Method 4: SMTTD (FFT)
        j=j+1
        opts=[];
        alpha=[1;1;1;1];
        alphasum = sum( alpha);
        alpha = alpha./ alphasum;
        opts.alpha=alpha;
        fall = [1e-4,1e-8,1e-8];
        opts.gamma =fall(i)*alpha;
        opts.beta = 10^(-4);
        opts.tau=[1.1,1];
        opts.max_iter=500;
        tolall = [1e-5,1e-3,1e-3];
        opts.tol = tolall(i);
        opts.T0=T;
        opts.Omega = Omega;
        opts.astv=1;
        opts.lambda=10;
        opts.weight=[0.5,0.5,0,0.5];
        transform.L = @fft; transform.l = 1; transform.inverseL = @ifft;
        tic;
        Xall{j} = SMTTD_LRTC4D(Y,opts,transform);
        Time(i,j)=toc;
        [PSNR(i,j),RSE(i,j),SSIM(i,j)]=quality_ll_order4(Xall{j}.*255,T.*255);
        %% Method 5: SMTTD (FFT-Data)
        j=j+1
        opts.gamma =10^(-4)*alpha;
        opts.beta = 10^(-4);
        tolall = [1e-3,2e-3,2e-3];
        opts.tol = tolall(i);
        [Xall{j},err] = SMTTD_LRTCU4D(Xall{j-1} ,Y,opts,transform);
        Time(i,j)=toc+Time(i,j-1);
        [PSNR(i,j),RSE(i,j),SSIM(i,j)]=quality_ll_order4(Xall{j}.*255,T.*255);
        %% Method 6: SMTTD (DCT)
        j=j+1
        opts=[];
        alpha=[1;1;1;1];
        alphasum = sum( alpha);
        alpha = alpha./ alphasum;
        opts.alpha=alpha;
        fall = [1e-4,1e-8,1e-8];
        opts.gamma =fall(i)*alpha;
        opts.beta = 10^(-3);
        opts.tau=[1.1,1];
        opts.max_iter=500;
        tolall = [1e-5,1e-3,1e-3];
        opts.tol = tolall(i);
        opts.T0=T;
        opts.Omega = Omega;
        opts.astv=1;
        opts.lambda= 10;
        opts.weight=[1,1,0,1];
        transform.L = @dct; transform.l = 1; transform.inverseL = @idct;
        tic;
        Xall{j} = SMTTD_LRTC4D(Y,opts,transform);
        Time(i,j)=toc;
        [PSNR(i,j),RSE(i,j),SSIM(i,j)]=quality_ll_order4(Xall{j}.*255,T.*255);
        %% Method 7: SMTTD (DCT-Data)
        j=j+1
        opts.gamma =10^(-4)*alpha;
        tolall = [1.5e-2,6e-3,2e-3]; %for foreman
        %         tolall = [8e-3,3e-3,2e-3]; %for suize
        opts.tol = tolall(i);
        [Xall{j},err] = SMTTD_LRTCU4D(Xall{j-1} ,Y,opts,transform);
        Time(i,j)=toc + Time(i,j-1);
        [PSNR(i,j),RSE(i,j),SSIM(i,j)]=quality_ll_order4(Xall{j}.*255,T.*255);
        
        savePath= [filename,'\', data_name,'_SR',num2str(SR),'j',num2str(j),'.mat'];
        save(savePath, 'PSNR','RSE','SSIM','Time','Xall','Y','Omega','SR','T');  %保存mat文件
    end
end

