clc;
clear all;
close all;
addpath(genpath('.\Test_data'));
addpath(genpath('.\mybib'));
filename = ['results_realdata_HSI'];
mkdir(filename);
%en=  {[370,75], [400,88], 6 , 2 }; 
en= {[386,73], [400,86], 10 , 2 };
for data_num =1:1:1
    switch data_num
        case 1
            load EO01_deadline.mat
            data_name = 'HSI_deadline';   
    end
    Y = Y_tensor0;
    Nway =size(Y);
    Omega = known;
    O = zeros(Nway);
    O(Omega) = 1;
    j=0;
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
%     opts.T0 = T; % original image
    % Choose transform
    transform.L = @fft; transform.l = 1; transform.inverseL = @ifft;
    % transform.L = @dct; transform.l = 1; transform.inverseL = @idct;
    % main loop
    tic;
    [Xall{j},err] = MTTD_LRTC(Y, Omega, opts, transform);
    % Computing the indexes
    Time(j)=toc;
%     [PSNR(j),RSE(j),SSIM(j)]=quality_ll(Xall{j}.*255,T.*255);
   %% Method 2： MTTD (FFT-Data)
    j=j+1;
    tic;
    [Xall{j},err] = MTTD_LRTCU(Xall{j-1}, Y, Omega, opts, transform);%
    % Computing the indexes
    Time(j)=toc+ Time(j-1);
%     [PSNR(j),RSE(j),SSIM(j)]=quality_ll(Xall{j}.*255,T.*255);
    %% Method 3： MTTD (DCT)
    j=j+1;
    % Choose transform
    transform.L = @dct; transform.l = 1; transform.inverseL = @idct;
    tic;
    [Xall{j},err] = MTTD_LRTC(Y, Omega, opts, transform);
    % Computing the indexes
    Time(j)=toc;
%     [PSNR(j),RSE(j),SSIM(j)]=quality_ll(Xall{j}.*255,T.*255);
    %% Method 4： MTTD (DCT-Data)
    j=j+1;
    tic;
    [Xall{j},err] = MTTD_LRTCU(Xall{j-1}, Y, Omega, opts, transform);%
    % Computing the indexes
    Time(j)=toc+ Time(j-1);
%     [PSNR(j),RSE(j),SSIM(j)]=quality_ll(Xall{j}.*255,T.*255);
%% Method 5： SMTTD (FFT)
    j=j+1;
    opts=[];
    alpha=[0.001,1,1];
    alphasum = sum( alpha);
    alpha = alpha./ alphasum;
    opts.alpha=alpha;
    opts.gamma =10^-4*alpha;
    opts.beta = 10^-2;
    opts.tau=[1.1,1];
    opts.max_iter=500;
    opts.tol=1e-3;
%     opts.T0=T;
    opts.Omega = Omega;
    opts.astv=1;
    opts.lambda=1;
    opts.weight=[1,1,1];
    transform.L = @fft; transform.l = 1; transform.inverseL = @ifft;
    tic;
    [Xall{j}] = SMTTD_LRTC(Y,opts,transform);
    Time(j)=toc;
%     [PSNR(j),RSE(j),SSIM(j)]=quality_ll(Xall{j}.*255,T.*255);
    %% Method 6： MTTD (FFT-Data)
    j=j+1;
    %     opts.tol=1e-3;
    opts.tol=5e-4;
    tic;
    [Xall{j},err] = SMTTD_LRTCU(Xall{j-1},Y,opts,transform);
    Time(j)=toc+ Time(j-1);
%     [PSNR(j),RSE(j),SSIM(j)]=quality_ll(Xall{j}.*255,T.*255);
    %% Method 7： SMTTD (DCT)
    j=j+1;
    opts=[];
    alpha=[0.001,1,1];
    alphasum = sum( alpha);
    alpha = alpha./ alphasum;
    opts.alpha=alpha;
    opts.gamma =10^-4*alpha;
    opts.beta = 10^-2;
    opts.tau=[1.1,1];
    opts.max_iter=500;
    opts.tol=1e-3;
%     opts.T0=T;
    opts.Omega = Omega;
    opts.astv=1;
    opts.lambda=1;
    opts.weight=[1,1,1];
    transform.L = @dct; transform.l = 1; transform.inverseL = @idct;
    tic;
    [Xall{j}] = SMTTD_LRTC(Y,opts,transform);
    Time(j)=toc;
%     [PSNR(j),RSE(j),SSIM(j)]=quality_ll(Xall{j}.*255,T.*255);
    %% Method 8： MTTD (DCT-Data)
    j=j+1;
    %     opts.tol=1e-3;
    opts.tol=5e-4;
    tic;
    [Xall{j},err] = SMTTD_LRTCU(Xall{j-1},Y,opts,transform);
    Time(j)=toc+ Time(j-1);
%     [PSNR(j),RSE(j),SSIM(j)]=quality_ll(Xall{j}.*255,T.*255);
    save([filename,'\',data_name,'.mat'],'Y','Xall','Omega');
    
    %% imshow the recovery image by different methods
    close all;
    mkdir([filename,'\',data_name]);
    figure
    num = [1:length(Xall)];
%     imshow(Y(:,:,2),'border','tight');
    enlarge(Y(:,:,1),[en{1}], [en{2}],en{3},en{4});
    saveas(gcf,[filename,'\',data_name,'\Y','.png']);
    for i=1:length(num)
        figure
%         imshow(Xall{num(i)}(:,:,2),'border','tight');
        enlarge(Xall{num(i)}(:,:,1),[en{1}], [en{2}],en{3},en{4});
        saveas(gca,[filename,'\',data_name,'\Xall',num2str(num(i)),'.png']);
    end
end
