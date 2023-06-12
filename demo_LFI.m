clc;
clear all;
close all;
addpath(genpath('.\Test_data'));
addpath(genpath('.\mybib'));
filename = ['results_LFI'];
mkdir(filename);
for data_num = 1:1:1
    switch data_num
        case 1
            load Flowers4D.mat
            data_name = 'Flowers4D';
        case 2
            data_name = 'Bikes4D';
            load Bikes4D.mat
    end
    Nway =size(T);
    N = numel(Nway);
    %% Generate Omega data
    for i = 1:2:3
        j=0;
        SRall = [0.01,0.05,0.1];
        SR = SRall(i);                    % Sample ratio (SR), e.g. 0.1 = 10% known samples
        mr = (1-SR)*100;                % Missing ratio (mr);
        Y=zeros(Nway);
        % Generate known data
        P = round(SR*prod(Nway));
        Omega = randsample(prod(Nway),P);
        Y(Omega)=T(Omega);
        
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
        %% main loop
        tic;
        [Xall{j}] = MTTD_LRTC(Y, Omega, opts, transform);
        %% 指标计算
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
        % 指标计算
        Time(i,j)=toc;
        [PSNR(i,j),RSE(i,j),SSIM(i,j)]=quality_ll_order4(Xall{j}.*255,T.*255);
        %% Method 3: MTTD (Data)
        j=j+1
        tic
        if PSNR(i,j-1)>=PSNR(i,j-2)
            U=Xall{j-1};
        else
            U=Xall{j-2};
            alpha=[1;1;1;1]; 
            alphasum = sum( alpha);
            alpha = alpha./ alphasum;
            opts.alpha = alpha;
            f = 1e-4;
            opts.gamma =f*alpha;
        end
        epsilonall = [1e-3,2e-3,2e-3];
        opts.epsilon = epsilonall(i);
        [Xall{j},err] = MTTD_LRTCU(U , Y, Omega, opts, transform);
        % 指标计算
        Time(i,j)=toc + Time(i,j-1);
        [PSNR(i,j),RSE(i,j),SSIM(i,j)]=quality_ll_order4(Xall{j}.*255,T.*255);
        %% Method 4: SMTTD (FFT)
        j=j+1
        opts=[];
        alpha=[1;1;1;1];% for SR=0.05
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
%         opts.lambda=0.1/SR;
        opts.weight=[0.1,0.1,1,1];
        transform.L = @fft; transform.l = 1; transform.inverseL = @ifft;
        tic;
        Xall{j} = SMTTD_LRTC4D(Y,opts,transform);
        Time(i,j)=toc;
        [PSNR(i,j),RSE(i,j),SSIM(i,j)]=quality_ll_order4(Xall{j}.*255,T.*255);
        %% Method 5: SMTTD (FFT-Data)
        j=j+1
        opts.gamma =10^(-4)*alpha;
        opts.beta = 10^(-4);
        tolall = [2e-4,2e-3,2e-3];
        opts.tol = tolall(i);
        [Xall{j},err] = SMTTD_LRTCU4D(Xall{j-1} ,Y,opts,transform);
        Time(i,j)=toc+ Time(i,j-1);
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
        tolall = [1e-5,1e-3,2e-3];
        opts.tol = tolall(i);
        opts.T0=T;
        opts.Omega = Omega;
        opts.astv=1;
        opts.lambda=10;
%         opts.lambda=0.1/SR;
        opts.weight=[1,1,1,1];
        transform.L = @dct; transform.l = 1; transform.inverseL = @idct;
        tic;
        Xall{j} = SMTTD_LRTC4D(Y,opts,transform);
        Time(i,j)=toc;
        [PSNR(i,j),RSE(i,j),SSIM(i,j)]=quality_ll_order4(Xall{j}.*255,T.*255);
        %% Method 6: SMTTD (DCT-Data)
        j=j+1
        opts.gamma =10^(-4)*alpha;
        tolall = [1.5e-2,3e-3,3e-3]; % Bikes
        %        tolall = [1e-3,3e-3,3e-3];
        opts.tol = tolall(i);
        [Xall{j},err] = SMTTD_LRTCU4D(Xall{j-1} ,Y,opts,transform);
        Time(i,j)=toc + Time(i,j-1);
        [PSNR(i,j),RSE(i,j),SSIM(i,j)]=quality_ll_order4(Xall{j}.*255,T.*255);

        filename1 = [filename, '\',data_name,'SR',num2str(SR)];
        mkdir(filename1);
        for n = 1 : length(Xall{j})
            X = Xall{n};
            save([filename1,'\Xa',num2str(j),'.mat'],'X');
        end 
        save( [filename1,'\','_j',num2str(j),'.mat'], 'PSNR','RSE','SSIM','Time','Y','Omega','SR','T');    
    end
end