function [X,E,err] = MTTD_RTPCA(T, opts, transform)
if ~exist('opts', 'var')
    opts = [];
end
if isfield(opts, 'T0');          T0= opts.T0;                 end
if isfield(opts, 'lambda');      lambda = opts.lambda;        end
if isfield(opts, 'alpha');       alpha= opts.alpha;           end
if isfield(opts, 'gamma');        gamma = opts.gamma;            end
if isfield(opts, 'max_gamma');       max_gamma = opts.max_gamma;          end
if isfield(opts, 'beta');       beta = opts.beta;          end
if isfield(opts, 'rho');         rho = opts.rho;              end
if isfield(opts, 'maxIter');     maxIter = opts.maxIter;      end
if isfield(opts, 'tol');     tol = opts.tol;      end

dim = size(T);
X=zeros(dim);
E=zeros(dim);
M = cell(ndims(T), 1);
for i=1:ndims(T)
    M{i}=zeros(dim);
    V{i}=zeros(dim);
end
Y=E;
Xtemp=M;
N = ndims(T);
if nargin < 3
    % fft is the default transform
    transform.L = @fft; transform.l = 1; transform.inverseL = @ifft;
end
% PSNR1=quality_ll(X.*255,T0.*255);
err =0;
chg =0;
for k = 1:maxIter
    gammasum = sum(gamma(:));
    tau = alpha./ gamma;
    if isfield(opts, 'T0')
        err(k) = RSE(T0(:),X(:));
        if mod(k,2) == 0
            fprintf('MTTD: iterations = %d   chg=%f RSE=%f\n', k, chg(k-1),err(k-1));
        end
    end
    Xsum = 0;
    Xold = X;
    Eold = E;
    %% Update M
    for i = 1:N
        Mtemp = Unfold_MTTD_RTPCA(X+V{i}/gamma(i), dim, i);
        M{i} = Fold_Square(prox_tnn(Mtemp, tau(i),transform), dim, i);
        Xtemp{i} = M{i}-V{i}/gamma(i);
        Xsum = Xsum+Xtemp{i}*gamma(i);    
    end
    
    %% Update X
    X = (Xsum + beta*(T+Y/beta-E))/(gammasum + beta);
    %% Update E
    E = prox_l1(T-X+Y/beta,lambda/beta);
     
    %% UpdateV
    for i = 1:N
       V{i}=V{i}-gamma(i)*(M{i}-X);
    end
     Y= Y+ beta*(T-E-X);
    gamma = min(rho(1)*gamma,max_gamma); 
    beta=rho(2)*beta;
    
    %% check the convergence
    dT = T-X-E;
    chg(k) =norm(X(:)-Xold(:))/norm(Xold(:));
    if  (k>=90)&&(chg(k) < tol)
        break;
    end 
end
