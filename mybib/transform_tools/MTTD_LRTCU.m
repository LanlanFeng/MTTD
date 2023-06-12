function [X,err] = MTTD_LRTCU(UU, T, Omega, opts, transform)
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
normT = norm(T(:));
dim = size(T);
M = cell(ndims(T), 1);
for i=1:ndims(T)
    M{i}=zeros(dim);
end
Y=M;
Xtemp=M;
N = ndims(T);
if nargin < 5
    % fft is the default transform
    transform.L = @fft; transform.l = 1; transform.inverseL = @ifft;
end
% transform.L = @dct; transform.l = 1; transform.inverseL = @idct
err=0;
for i=1:ndims(T)
    Utemp = Unfold_Square(UU, dim, i);
    O = tenmat(Utemp,[3]); % unfolding
    O = O.data;
    %     if size(O,1)<=size(O,2)
    [U0,D0,V0]= svd(O,'econ');
    %     else
    %         U0=zeros(size(O,1),size(O,1));
    %         [U0(:,1:size(O,2)),D0,V0]= svd(O,'econ');
    %     end
    U{i}=U0;
end

% PSNR1=quality_ll(X.*255,T0.*255);
for k = 1:maxIter
    gammasum = sum(gamma);
    tau = alpha./ gamma;
    if isfield(opts, 'T0')
        %         PSNR1=quality_ll(X.*255,T0.*255);
        err(k) = RSE(T0(:),X(:));
        if mod(k,2) == 0
            %        fprintf('MTTD-Data: iterations = %d   difference=%f PSNR=%f\n', k, errList(k-1),PSNR1);
            fprintf('MTTD-Data: iterations = %d   difference=%f RSE=%f\n', k, errList(k-1),err(k-1));
        end
    end
    
    Msum = 0;
    Ysum = 0;
    Xsum = 0;
    %% Update M
    for i = 1:N
        Mtemp = Unfold_Square(X+Y{i}/gamma(i), dim, i);
        n3 = size(Mtemp,3);
        L=U{i};transform.l = 1; transform.L = L;
        M{i} = Fold_Square(prox_utnn(U{i},Mtemp, tau(i)), dim, i);
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
    if (k>=50)&&(errList(k) < epsilon) &&(errList(k)>1e-8)%may change accroding to different data
%     if (k>=50)&&( err(k)<0.9)&&(errList(k) < epsilon) &&(errList(k)>1e-8)%may change accroding to different data
        errList = errList(1:k);
        break;
    end
%     if (ndims(T)==3)&&(dim(3) ==3)
%         imshow(X)
%     end
end