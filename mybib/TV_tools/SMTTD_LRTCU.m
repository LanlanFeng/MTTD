function [X,err] = SMTTD_LRTCU(UU,T,opts, transform)
if ~exist('opts', 'var')
    opts = [];
end
if isfield(opts, 'tol');         tol = opts.tol;              end
if isfield(opts, 'max_iter');    max_iter = opts.max_iter;    end
if isfield(opts, 'lambda');      lambda = opts.lambda;        end
if isfield(opts, 'alpha');       alpha = opts.alpha;          end
if isfield(opts, 'gamma');        gamma = opts.gamma;          end
if isfield(opts, 'beta');        beta = opts.beta;          end
if isfield(opts, 'tau');         tau = opts.tau;              end
if isfield(opts, 'weight')        weight = opts.weight;       end
if isfield(opts, 'Omega')        Omega = opts.Omega;       end
if isfield(opts, 'T0');          T0= opts.T0;                 end
%% Initialization X,M,L;
dim = size(T);
N = length(dim);
X = T;
M = cell(N,1);
for i=1:N
    M{i} = zeros(dim); %% the auxiliary tensor
end
R  = prod(dim);
L = zeros(3*R, 1); %% the auxiliary variable
V=M;
Y=L;

for i=1:ndims(T)
    Utemp = Unfold_Square(UU, dim, i);
    O = tenmat(Utemp,[3]); % unfolding
    O = O.data;
    %     if size(O,1)<=size(O,2)
    [U0,D0,V0]= svd(O,'econ');
    %     else
    %     U0=zeros(size(O,1),size(O,1));
    %     [U0(:,1:size(O,2)),D0,V0]= svd(O,'econ');
    %     end
    % %     [U0,D0,V0]= svd(O);
    U{i}=U0;
end
%% Initialization: Compute D^* D
h    = dim(1);
w    = dim(2);
d    = dim(3);
Eny_x   = ( abs(psf2otf([-1; +1], [h,w,d])) ).^2  ;
Eny_y   = ( abs(psf2otf([-1, +1], [h,w,d])) ).^2  ;
Eny_z   = ( abs(psf2otf([+1, -1], [w,d,h])) ).^2  ;
Eny_z   =  permute(Eny_z, [3, 1 2]);
% denom1  =  Eny_x + Eny_y ;
if isfield(opts, 'astv')
    denom1  =  (weight(1)^2)*Eny_x +  (weight(2)^2)*Eny_y +  (weight(3)^2)*Eny_z;
else
    denom1  =  weight(1)*Eny_x +  weight(2)*Eny_y +  weight(3)*Eny_z;
end
%% Main Loop
Out.Res=[];Out.RSE=[];Out.PSNR=[];
iter = 0;
err =0;
for iter = 1 : max_iter
    Xold =X;
    
    %% Updating auxiliary variable M about X
    Msum =0;
    Vsum =0;
    for  i = 1:N
        %      transform.L=U{i};transform.l = 1; transform.L = L;
        M{i} = Fold_Square(prox_utnn(U{i},Unfold_Square(X-V{i}/gamma(i), dim, i), alpha(i)/gamma(i)), dim, i);
        Msum = Msum + gamma(i)*M{i};
        Vsum = Vsum + V{i};
    end
    
    %% Updating low-rank component X
    temp_X =Msum+Vsum;
    if isfield(opts, 'astv')
        diffT_L = diffT3( beta*L+Y, dim ,weight);%AS
    else
        diffT_L = diffT3_2( beta*L+Y, dim ,weight);%IS
    end
    numer1   = reshape( diffT_L + temp_X(:), dim);
    gammasum = sum(gamma);
    X = real( ifftn( fftn(numer1) ./ (beta*denom1 + gammasum) ) );
    X(Omega) = T(Omega);
    
    %% Updating auxiliary variable L
    if isfield(opts, 'astv')
        diff_x = diff3(X, dim, weight);
        L      = softThres( diff_x - Y/beta, lambda/beta); %AS
    else
        diff_x = diff3_2(X, dim, weight);
        L      = softThres21( diff_x - Y/beta, lambda/beta, R); %IS
    end
    %% Updating dual variable V and Y
    for  i = 1:N
        V{i} = V{i}+ gamma(i)*(M{i}-X);
    end
    Y= Y+ beta*(L-diff_x);
    gamma=tau(1).*gamma;
    beta=tau(2)*beta;
    %% check the convergence
    
    res=norm(X(:)-Xold(:))/norm(Xold(:));
    Out.Res = [Out.Res,res];
    %      if iter==1 || mod(iter, 2) == 0
    if isfield(opts, 'T0')
        err(iter) = RSE(T0(:),X(:));
        %         PSNR1=quality_ll(X.*255,XT.*255);
        %     Out.PSNR = [Out.PSNR,PSNR1];
        if mod(iter,2) == 0
            %        fprintf('SMTTD-Data: iterations = %d   difference=%f PSNR=%f\n', k, errList(k-1),PSNR1);
            fprintf('SMTTD-Data: iterations = %d   difference=%f RSE=%f\n', iter, res,err(iter-1));
        end
    end
        if (iter>=40)&&(res < tol)%may change accroding to different data
%     if (err(iter)<0.9)&&(res < tol)
        break;
    end
    %     iter
%     if (ndims(T)==3)&&(dim(3) ==3)
%         imshow(X)
%     end
    
end
