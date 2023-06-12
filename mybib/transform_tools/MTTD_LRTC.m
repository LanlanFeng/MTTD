function [X,err,errList] = MTTD_LRTC(T, Omega, opts, transform)

if ~exist('opts', 'var')
    opts = [];
end
if isfield(opts, 'T0');          T0= opts.T0;                 end
if isfield(opts, 'alpha');       alpha= opts.alpha;           end
% if isfield(opts, 'beta2');        beta2 = opts.beta2;            end
if isfield(opts, 'gamma');       gamma = opts.gamma;          end
if isfield(opts, 'rho');         rho = opts.rho;              end
if isfield(opts, 'maxIter');     maxIter = opts.maxIter;      end
if isfield(opts, 'epsilon');     epsilon = opts.epsilon;      end

X = T;
dim = size(T);
M = cell(ndims(T), 1);
for i=1:ndims(T)
    M{i}=zeros(dim);
end
Y=M;
Xtemp=M;
N = ndims(T);
if nargin < 4
    % fft is the default transform
    transform.L = @fft; transform.l = 1; transform.inverseL = @ifft;
end
err =0;
for k = 1:maxIter
    gammasum = sum(gamma);
    tau = alpha./ gamma;
    if isfield(opts, 'T0')
        %       PSNR1 = quality_ll(X(:).*255,T0(:).*255);
        err(k) = RSE(T0(:),X(:));
        if mod(k,2) == 0
            %        fprintf('MTTD: iterations = %d   difference=%f PSNR=%f\n', k, errList(k-1),PSNR1);
            fprintf('MTTD: iterations = %d   difference=%f RSE=%f\n', k, errList(k-1),err(k-1));
        end
    end
    Msum = 0;
    Ysum = 0;
    Xsum = 0;
    %% Update M
    for i = 1:N
        Mtemp = Unfold_Square(X+Y{i}/gamma(i), dim, i);
        M{i} = Fold_Square(prox_tnn(Mtemp, tau(i),transform), dim, i);
        Xtemp{i} = M{i}-Y{i}/gamma(i);
        Xsum = Xsum+Xtemp{i}*gamma(i);
    end
    Xlast = X;
    %% Update X
    X =Xsum/gammasum;
    X(Omega) = T(Omega);
    %% Update Y
    for i = 1:N
        Y{i}=Y{i}-gamma(i)*(M{i}-X);
    end
    gamma=gamma*rho;
    errList(k) = norm(X(:)-Xlast(:)) / norm(Xlast(:));
%     if (k>=45)&&( err(k)<0.9)&&(errList(k) < epsilon) % may change accroding to different data
    if (k>=45)&&(errList(k) < epsilon) % may change accroding to different data
        %     if (k>240)
        errList = errList(1:k);
        break;
    end
    
%     if (ndims(X)==3)&&(dim(3) ==3)
%         imshow(X)
%     end
    %     if (ndims(X)==4)&&(dim(3) ==13)
    %         imshow(X(:,:,[4,3,2],3));
    %     end
    %     if (ndims(X)==2)&&(dim(3) ==3)
    %         imshow(X(:,:,:,2));
    %     end
end