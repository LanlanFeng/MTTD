clc;
clear all;
close all;
addpath(genpath('.\Test_data'));
addpath(genpath('.\mybib'));
filename = ['results_colorimage_imshow'];
mkdir(filename);
for data_num = 1:1:3
    switch data_num
        case 1
            load('peppers.mat');
            data_name = 'peppers0.01';
            SR=0.01;
        case 2
            load('lena.mat');
            data_name = 'lena0.05';
            SR=0.05;
        case 3
            load('house.mat');
            data_name = 'house0.1';
            SR=0.1;
        case 4
            load('inpainting1.mat')
            data_name = 'barbara_word';
        case 5
            load('inpainting2.mat')
            data_name = 'lena_gridlines';
        case 6
            load('sails_deadline.mat')%tulips_block
            data_name = 'sails_deadline';
    end
    
    %% Generate Omega data
    % SR =0.05;%%change the sample ratio
    if data_num<=3
        Nway =size(T);
        N = numel(Nway);
        Y=zeros(Nway);
        P = round(SR*prod(Nway));
        Omega = randsample(prod(Nway),P);
        Y(Omega)=T(Omega);
    else
        Y = Y_tensor0;
        T = Y_tensorT;
        Omega = known;
        Nway =size(T);
        N = numel(Nway);
    end
    %     O = zeros(Nway);
    %     O(Omega) = 1;
    
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
    opts.lambda=10;
    opts.weight=[1,1,0];
    transform.L = @dct; transform.l = 1; transform.inverseL = @idct;
    tic;
    [Xall{j}] = SMTTD_LRTC(Y,opts,transform);
    Time(j)=toc;
    [PSNR(j),RSE(j),SSIM(j)]=quality_ll(Xall{j}.*255,T.*255);
    %% Method 5： MTTD (DCT-Data)
    j=j+1;
    opts.weight=[0.5,0.5,0];
%     opts.weight=[1,1,0];
    %         opts.tol=2e-3;
    tic;
    [Xall{j},err] = SMTTD_LRTCU(Xall{j-1},Y,opts,transform);
    Time(j)=toc+ Time(j-1);
    [PSNR(j),RSE(j),SSIM(j)]=quality_ll(Xall{j}.*255,T.*255);
    
    %% save  mat file
    save([filename,'\',data_name,'.mat'],'T','Y','Xall','Omega');
    
    %% imshow the recovery image by different methods
    data_numall = 6; % # the samples
    num = [1:3,5];
    colnum = length(num)+1; % # the columns
    FS=15; % FontSize
    label = ['bcdef'];
    
    subplot = @(m,n,p) subtightplot (m,n,p,[0.0028 0.0028], [0.02 0.015], [0.002 0.002]);% the parameter of the suptitle
    %First []: gaps between subgraphs
    %Second[]: blank space for the bottom/top of total figure
    %Third []: blank space for the left/right of  the total figure
    subplot(data_numall,colnum,(colnum)*(data_num-1)+1)
    imshow(Y);
    if (data_num==data_numall)
        xlabel('(a)','VerticalAlignment','baseline','FontSize',FS,'FontName','Times');
    end
    title([data_name],'FontSize',FS,'FontName','Times')
    for i=1:1:length(num)
        subplot(data_numall,colnum,(colnum)*(data_num-1)+i+1)
        imshow(Xall{num(i)});
        a=quality_ll(Xall{num(i)}.*255,T.*255);
        a= sprintf('%6.4f',a) ;
        title([num2str(a),'dB'],'FontSize',FS,'FontName','Times')
        if (data_num==data_numall)
            xlabel(['(',label(i),')'],'VerticalAlignment','baseline','FontName','Times','FontSize',FS);
        end
    end
    %     %% save the pdf file
    %     set(gcf,'Units','Inches');
    %     pos = get(gcf,'Position');
    %     set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    %     filename = 'colorimageimpainting';
    %     print(gcf,filename,'-dpdf','-r0')
    %     close all;
    %     mkdir([filename,'\',data_name]);
    %     figure
    %     num = [1:5];
    %     imshow(Y,'border','tight');
    %     saveas(gcf,[filename,'\',data_name,'\Y','.png']);
    %     for i=1:length(num)
    %         figure
    %         imshow(Xall{num(i)},'border','tight');
    %         saveas(gca,[filename,'\',data_name,'\Xall',num2str(num(i)),'.png']);
    %     end
end

