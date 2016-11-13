function VecFld = MVFC_test(varargin)
% Test version for the  L * L neighbours and the width of the Gaussian filter H 
% VecFld = MVFC_test(Vx,Vy,sigma0,L,H)

% MVFC: Modified Vector Field Correction method for PIV outlier detection
% Syntax:
%   VecFld = MVFC(Vx,Vy)
%   VecFld = MVFC(Vx,Vy,Vz)
%   VecFld = MVFC(Vx,Vy,sigma0)
% Input:
%   Vx,Vy : The velocity field in gridded.
%   sigma0 : A small value (coefficient) to void the sigma close to zero, and can be adjusted in [0,0.01] as a threshold. Default value is 0.001.
% Output:
%   VecFld: A structure type value which contains X, V, beta, V, P,VFCIndex.
%       Where V denotes the transformation value of X by VFC, i.e the
%       processed velocity field of U. P is the posterior probability and VFCIndex is the
%       indexes of inliers which found by VFC.
%
%   See also:: smoothn(), which must be required in VFC.
%   Refs:
%        1. Garcia, D. (2010). Robust smoothing of gridded data in one and higher dimensions with missing values. Computational statistics & data analysis, 54(4), 1167-1178.
%        2. Garcia, D. (2011). A fast all-in-one method for automated post-processing of PIV data. Experiments in fluids, 50(5), 1247-1259.
%        3. Ma, J., Zhao, J., Tian, J., Yuille, A. L., & Tu, Z. (2014). Robust point matching via vector field consensus. IEEE Transactions on Image Processing, 23(4), 1706-1721.
%        4. Lee, Y., Yang, H., & Yin, Z. (2016). A robust vector field correction method via a mixture statistical model of PIV signal. Experiments in Fluids, 57(3), 1-20.
%        5. Lee, Y., Yang, H. (2016). Statistical outlier detection for particle image velocimetry data using a local noise variance estimate. Submitting to Measurement Science and Technology .
%
% % Example
% % Cellular vortical flow
%   [x,y] = meshgrid(linspace(0,1,24));
%   Vx = cos(2*pi*x+pi/2).*cos(2*pi*y);
%   Vy = sin(2*pi*x+pi/2).*sin(2*pi*y);
%   Vx = Vx + (0.01)*randn(24,24); % adding Gaussian noise
%   Vy = Vy + (0.01)*randn(24,24); % adding Gaussian noise
%   I = randperm(numel(Vx));
%   Vx(I(1:60)) = (rand(60,1)-0.5)*5; % adding outliers
%   Vy(I(1:60)) = (rand(60,1)-0.5)*5; % adding outliers
%   Vx(I(61:100)) = NaN; % missing values
%   Vy(I(61:100)) = NaN; % missing values
%   VecFld = VFC(Vx,Vy,0.00);  %VFC Method
%   Vx_VFC = VecFld.V(:,:,1);Vy_VFC = VecFld.V(:,:,2);VFC_Flag = VecFld.VFCIndex;% Output
%   subplot(121), quiver(x,y,Vx,Vy,2.5), axis square
%   title('Noisy velocity field')
%   subplot(122), quiver(x,y,Vx_VFC,Vy_VFC), axis square
%   title('Processed velocity field')

%%  Check the input
% number of argument input:  nargin
% variable of the arguments: varargin{k}
if nargin<1;demo();return;end;% run the default demo example.

Compont1 = varargin{1};   % Get the first component
[Nx,Ny,~]=size(Compont1); % Get the 2-D size
for k = 1: nargin
  if size(varargin{k}) == size(Compont1) % Add another component
      U(:,k)= varargin{k}(:);
%       U(:,k) = U(:,k) + max(abs(U(:,k)))*(0.001)*randn(size(U(:,k))); % Add a small noise
  else 
      break;
  end
  C = k; % Recode the Component Number as C
end

% Get other settings 
DefaltSigma0Coef = 0.001; % a default variable constant
if  C == nargin 
        sigma0Coef = DefaltSigma0Coef;
else if    isscalar(varargin{k})
    sigma0Coef = varargin{k};
    else 
      error('MATLAB:MVFC:NotScaleCoefficient',...
            'Please check the input.')
    end
end

L = varargin{k+1};
H = varargin{k+2};
%- Set the default minimum value restrition on sigma0: 
if sigma0Coef<0 || sigma0Coef>=1
    sigma0Coef = DefaltSigma0Coef;
     warning('MATLAB:MVFC:DefaultPara',...
            'unproper sigma0 adopted')
end

% Change the data into 3-D Matrix form
sizX = size(Compont1); 
U = reshape(U,[sizX, C]);

%% begin the code
% fprintf('Starting Modified Vector Field Correction:\n');

% Initialization
V = smoothvfc(U,sizX,C);         % estimate the initial reference-field with a large smooth kernel
I = isfinite(U); U(~I) = V(~I);  % set a random value to replace missing value
IndMiss = true(sizX); for i=1:C;  IndMiss = IndMiss & ~I(:,:,i);end
ConstSigma2(C)=0; sigma0 = 0*U; Sigma2 = 0*U;  % initial global variance on each component
                                               % A small sigma0, avoid the
                                               % Sigma2 close to zere in a
                                               % iteration

for i=1:C % For each component, estimate the variance and set the small bias sigma0
ConstSigma2(i)=sum(sum((U(:,:,i)-V(:,:,i)).^2))/length(find(I(:,:,i))); % The globle variance estimation
Sigma2(:,:,i) = ConstSigma2(i);
sigma0(:,:,i) = ConstSigma2(i)*sigma0Coef;
% ConstSigma2(i)*sigma0Coef % For debug display
% Temp=abs(U(:,:,i)-V(:,:,i));a(i) = 2*max(Temp(:));
end
a=1.8*sqrt(ConstSigma2); % Set the outlier distribute range

%- Outlier ratio, threshold value, maxmumum iteration number and error 
gamma=0.90;tau=0.25; MaxIter=21; ecr=2*10^-8; iter=1;  tecr=1; 

%% Iteration of EM algorithm
NewP=zeros(Nx,Ny); OutlierNum = round((1-gamma)*(Nx*Ny)); OutlierRatio =1- gamma;
while (iter<MaxIter) && (tecr > ecr) && (max(Sigma2(:)) > 1e-10)
    % E-step: Update P and calculate the changes (for termination)
    P=getPosterior(U,V, Sigma2, gamma, a, I, sigma0);
    tecr = max(max(abs(P-NewP))); NewP = P;
    fprintf('iterate: %dth,  Outlier Number (Ratio): %d (%f), Missing Number: %d, maxSigma2=%f\n', iter, OutlierNum, OutlierRatio,sum(IndMiss(:)), max(Sigma2(:)));
    
    % M-step.
    %- Update V on each component
    warning off all
    for k = 1:size(V,3)
        Temp = U(:,:,k); TempMean = mean(Temp(:));
        [Z,S] = smoothn(U(:,:,k)-TempMean,P);      % To satisfy the conditon 
        if S<0.01 %The data should 
             [Z,S] = smoothn(U(:,:,k)-TempMean,P,0.01);
        end
      V(:,:,k)=Z+TempMean;
    end
    warning on all

    %- Update gamma
    numcorr = length(find((P(:) > tau)&~IndMiss(:))); OutlierNum = Nx*Ny-numcorr;
    gamma = C*numcorr/(sum(I(:))); OutlierRatio =1-gamma;
    if gamma > 0.95, gamma = 0.95; end
    if gamma < 0.05, gamma = 0.05; end
    
    %- Update Sigma2 
    for i =1:C
    SigmaTemp = ((P+0.00).*(U(:,:,i)-V(:,:,i)).^2);
%     sum(SigmaTemp(:))./sum(P(:))% For debug display
    h = fspecial('gaussian', [L,L], H); % Gaussian Kernel
%         h = fspecial('average', [16,16]);% Uniform Kernel
%         h = 1;                           % Without Kernel
    Sigma2(:,:,i) = imfilter(SigmaTemp,h,'replicate')./imfilter(P,h,'replicate');
    end
    
    %- iteration number updata
    iter=iter+1;
end

% Do a small smoothness and interpolation at last
V(:,:,1)=smoothn(V(:,:,1),(P > tau)&~IndMiss,0.05);V(:,:,2)=smoothn(V(:,:,2),(P > tau)&~IndMiss,0.05);

%% Output
VecFld.Y = U;
VecFld.V = V;
VecFld.P = P;
VecFld.Index = (P > tau);
disp('Removing outliers succesfully completed.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = getPosterior(U,V,Sigma2,gamma,a,I,Sigma0)
% getPosterior estimates the posterior probability
[sx,sy,d] = size(U); 
P = ones(sx,sy);
X = U-V; X(~I)=0;% The residual part

Sigma2 = Sigma2 + Sigma0; % Add a small value avoiding the zero dominator
P_inlier = exp(-X.^2./(2*Sigma2))./sqrt(2*pi*Sigma2);
for i = 1:d
    TempP = P_inlier(:,:,i)./(P_inlier(:,:,i)+(1-gamma)/(a(i)*gamma));
    P = P.*TempP;
end
%     P = (P-min(P(:)))/(max(P(:))-min(P(:)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function V=smoothvfc(Y,sizX,C)
for i=1:C
Temp = Y(:,:,i); Temp = Temp(:);Temp(find(isnan(Temp)))=[];
TempMean = mean(Temp(:)); % Calculate the mean of the valid component data

V(:,:,i)=smoothn(Y(:,:,i)-TempMean,1000)+TempMean;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function demo()
N = 32;
[x,y] = meshgrid(linspace(0,1,N));
Vx = cos(2*pi*x+pi/2).*cos(2*pi*y);
Vy = sin(2*pi*x+pi/2).*sin(2*pi*y);

rng('default');
rng(10);Vx = Vx + (0.005)*randn(N,N); % adding Gaussian noise
rng(20);Vy = Vy + (0.005)*randn(N,N); % adding Gaussian noise
rng(30);I = randperm(numel(Vx));
Vx(I(1:60)) = (rand(60,1)-0.5)*4; % adding outliers
Vy(I(1:60)) = (rand(60,1)-0.5)*4; % adding outliers
Vx(10,10) = NaN;
Vx(I(61:100)) = NaN; % missing values
Vy(I(61:100)) = NaN; % missing values

tic
VecFld = MVFC_test(Vx,Vy,0.001,16,10);% MVFC Method
toc

Vx_MVFC = VecFld.V(:,:,1);Vy_MVFC = VecFld.V(:,:,2);VFC_Flag = VecFld.Index;% Output
subplot(131), quiver(x,y,Vx,Vy,2.5), axis square;
xlim([-0.05,1.05]);ylim([-0.05,1.05]);box on;set(gca,'ytick',[]);set(gca,'xtick',[]);
title('Noisy velocity field')
subplot(132), quiver(x,y,Vx_MVFC,Vy_MVFC), axis square;
xlim([-0.05,1.05]);ylim([-0.05,1.05]);box on;set(gca,'ytick',[]);set(gca,'xtick',[]);
title('Processed velocity field by MVFC')
subplot(133), spy(~flipud(VFC_Flag)), axis square;
xlim([-0.05,N+0.05]);ylim([-0.05,N+0.05]);box on;set(gca,'ytick',[]);set(gca,'xtick',[]);
title('The outlier distribution by MVFC')
end
