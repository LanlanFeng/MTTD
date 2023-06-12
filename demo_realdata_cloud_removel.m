clc;
clear all;
close all;
addpath(genpath('.\Test_data'));
addpath(genpath('.\mybib'));
filename = ['results_cloudremovel'];
mkdir(filename);
for datanum =1:1:1
    switch datanum
        case 1
            load('cloudy_data.mat')%
            data_name = 'cloudy';
    end
    Y = Y_tensor0;
    %     T = Y_tensorT;
    Nway =size(Y);
    N = numel(Nway);
    Omega = known;
    %     O = zeros(Nway);
    %     O(Omega) = 1;
    
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
    opts.epsilon = 4e-3;
    %     opts.T0 = T; % original image
    % Choose transform
    transform.L = @fft; transform.l = 1; transform.inverseL = @ifft;
    %% main loop
    [Xall{j}] = MTTD_LRTC(Y, Omega, opts, transform);
    %% Method 2: MTTD (FFT-Data)
    j=j+1
    tic
    U=Xall{j-1};
    opts.gamma =1e-3*alpha;
    opts.epsilon = 2e-3;
    [Xall{j}] = MTTD_LRTCU(U , Y, Omega, opts, transform);
    %% Method 3: MTTD (DCT)
    j=j+1
    alpha=[1;1;1;1];
    alphasum = sum( alpha);
    alpha = alpha./ alphasum;
    opts.alpha = alpha;
    f = 1e-4;
    opts.gamma =f*alpha;
    opts.epsilon = 1e-3;
    transform.L = @dct; transform.l = 1; transform.inverseL = @idct;
    [Xall{j}] = MTTD_LRTC(Y, Omega, opts, transform);
    %% Method 4: MTTD (DCT-Data)
    j=j+1
    tic
    U=Xall{j-1};
    opts.gamma =1e-3*alpha;
    opts.epsilon = 2e-3;
    [Xall{j}] = MTTD_LRTCU(U , Y, Omega, opts, transform);
    %% Method 5: SMTTD (FFT)
    j=j+1
    opts=[];
    alpha=[1;1;1;1];
    alphasum = sum( alpha);
    alpha = alpha./ alphasum;
    opts.alpha=alpha;
    opts.gamma =1e-4*alpha;
    opts.beta = 10^(-3);
    opts.tau=[1.1,1];
    opts.max_iter=500;
    opts.tol = 4e-3;
    % opts.T0=T;
    opts.Omega = Omega;
    opts.astv=1;
    wlall=[10;50;100;500];
    opts.lambda=10;
    opts.weight=[1,1,0,1];
    transform.L = @fft; transform.l = 1; transform.inverseL = @ifft;
    Xall{j} = SMTTD_LRTC4D(Y,opts,transform);
    %% Method 6: SMTTD (FFT-Data)
    j=j+1
    opts.gamma =10^(-4)*alpha;
    opts.beta = 10^(-3);
    opts.tol = 2e-3;
    [Xall{j}] = SMTTD_LRTCU4D(Xall{j-1} ,Y,opts,transform);
    %% Method 7: SMTTD (DCT)
    j=j+1
    opts=[];
    alpha=[1;1;1;1];
    alphasum = sum( alpha);
    alpha = alpha./ alphasum;
    opts.alpha=alpha;
    fall = [1e-4,1e-8,1e-8];
    opts.gamma =fall(1)*alpha;
    opts.beta = 10^(-3);
    opts.tau=[1.1,1];
    opts.max_iter=500;
    opts.tol = 1e-3;
    % opts.T0=T;
    opts.Omega = Omega;
    opts.astv=1;
    opts.lambda=10;
    opts.weight=[1,1,0,1];
    transform.L = @dct; transform.l = 1; transform.inverseL = @idct;
    Xall{j} = SMTTD_LRTC4D(Y,opts,transform);
    %% Method 8: SMTTD (DCT-Data)
    j=j+1
    opts.gamma =10^(-4)*alpha;
    opts.tol = 1e-3;
    [Xall{j}] = SMTTD_LRTCU4D(Xall{j-1} ,Y,opts,transform);
    savePath= [filename,'\', data_name,'.mat'];
    save(savePath,'Xall','Y','Omega');
end
figure
label = ['bcdefghijklmno'];
data_numall = 1; % # the samples
FS=12; % FontSize
% en= {[10,180], [50,220], 3.5 , 2 };
% en= {[50,100], [100,160], 3 , 2 }; %frame 2
% en= {[50,130], [70,160], 5 , 2 }; %frame 2
en= {[50,130], [70,155], 6 , 2 }; %frame 2
num=[1:length(Xall)]; %the order of the method
colnum = length(num)+1;
data_num =1;
frame(data_num)=2;
subplot = @(m,n,p) subtightplot (m,n,p,[0.0028 0.0028], [0.16 0.001], [0.002 0.002]);
%First []: gaps between subgraphs
%Second[]: blank space for the bottom/top of total figure
%Third []: blank space for the left/right of  the total figure
subplot(1,colnum,colnum*(data_num-1)+1)
enlarge(Y(:,:,:,frame(data_num)),[en{1}], [en{2}],en{3},en{4});
if (data_num==data_numall)
    xlabel('(a)','VerticalAlignment','baseline','FontSize',FS,'FontName','Times');
end
for n=1:length(num)
    subplot(1,colnum,colnum*(data_num-1)+n+1)
    enlarge(Xall{n}(:,:,:,frame(data_num)),[en{1}], [en{2}],en{3},en{4});
    if (data_num==data_numall)
        xlabel(['(',label(n),')'],'VerticalAlignment','baseline','FontName','Times','FontSize',FS);
    end
end
