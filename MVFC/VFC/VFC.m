function VecFld = VFC(Vx,Vy,sigma0)
% VFC  Vector Field Correction method implementation
%   VECFLD = VFC(Vx,Vy,sigma0)
%   learning vector field from grided samples with outliers.
%   
% Input:
%   Vx, Vy : The velocity field in gridded.
%   sigma0 : A small value to void the sigma too close to zero, and can be adjusted in [0,1] as a threshold. Default value is 1.0.
%
%
% Output:
%   VecFld: A structure type value which contains X, V, beta, V, P,VFCIndex.
%       Where V denotes the transformation value of X by VFC, i.e the
%       processed velocity field of U. P is the posterior probability and VFCIndex is the
%       indexes of inliers which found by VFC.
%
%   See also:: smoothn(), which is required in VFC.
%   Ref1: Yong Lee (2015) A Robust Vector Field Correction Method via a Mixture Statistical Model of PIV signal.submited to Exp Fluids
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
if nargin<1;demo();return;end;% run the default example.

if nargin < 3
    sigma0 = 1.0;% default
end

if ~isequal(size(Vx),size(Vy))
        error('MATLAB:VFC:SizeMismatch',...
            'Arrays for Vx and Vy must have same size.')
elseif  (~isscalar(sigma0) || sigma0<0)
    error('MATLAB:VFC:IncorrectParameter',...
        'The sigma0 parameter must be a scalar >=0')
end

if sigma0<0 || sigma0>1
    sigma0 = 1.0;
end

sigma0 = sigma0/20;  

%% Initialization
gamma=0.90; lambda=1; theta=0.5; a=[]; MaxIter=20; ecr=2*10^-8;
[Nx,Ny]=size(Vx);D = 2;

%% begin the code
fprintf('Starting Vector Field Correction Learning:\n');
U(:,:,1) = Vx; U(:,:,2) = Vy;
% Initialization
V=smoothvfc(U); iter=1;  tecr=1; E=1;K=[1 1 1 ;1 -8 1;1 1 1]; % C=zeros(N,D);
I  =  isfinite(U); % missing value support
U(~I) = V(~I);       % give a random value 0
sigma2=sum(sum(sum((U(I)-V(I)).^2)))/(sum(I(:)));a=2.5*sigma2;

%% Iteration of EM algorithm
newE = [];newE2 = []; NewP=zeros(Nx,Ny); OutlierNum = round((1-gamma)*(Nx*Ny)); OutlierRatio =1- gamma;
while (iter<MaxIter) && (tecr > ecr) && (sigma2 > 1e-10) 
    % E-step. 
%     E_old=E;
    [P, E]=get_P(U,V, sigma2 ,gamma, a,I,sigma0);
%     E=E+lambda/2*sum(sum((conv2(V(:,:,1),K,'valid').^2+conv2(V(:,:,2),K,'valid').^2)));
%     newE = [newE,E];%   
%     tecr=abs((E-E_old)/E);
%     newE2 = [newE2,tecr];%
    tecr  = max (max(abs(P-NewP)));NewP=P;
    fprintf('iterate: %dth,  Outlier Number (Ratio): %d (%f), sigma2=%f\n', iter, OutlierNum, OutlierRatio, sigma2);

    % M-step. 
    % --Update V and sigma^2
    warning off all
    [Vtemp,lambda]=smoothn(complex(U(:,:,1),U(:,:,2)),P,'Initial',complex(V(:,:,1),V(:,:,2)));
    V(:,:,1)=real(Vtemp);V(:,:,2)=imag(Vtemp);
    warning on all

    Sp=sum(P(:));
    sigma2=sum(sum(P.*sum((U-V).^2, 3)))/(Sp*D);
    % --Update gamma
    numcorr = length(find(P(:) > theta)); OutlierNum = Nx*Ny-numcorr;
    gamma = 2*numcorr/(sum(I(:))); OutlierRatio =1-gamma;
    if gamma > 0.95, gamma = 0.95; end
    if gamma < 0.05, gamma = 0.05; end
    
    iter=iter+1;
end
% Do a small smoothness and interpolation at last
V(:,:,1)=smoothn(V(:,:,1),P,0.05);V(:,:,2)=smoothn(V(:,:,2),P,0.05);

%% Output
VecFld.Y = U;
VecFld.V = V;
VecFld.P = P;
VecFld.VFCIndex = (P > theta);
disp('Removing outliers succesfully completed.');
end

%%%%%%%%%%%%%%%%%%%%%%%%
function [P, E]=get_P(Y,V, sigma2 ,gamma, a,I,sigma0)
% GET_P estimates the posterior probability and part of the energy.
Y(~I)=0;V(~I)=0;
[Nx,Ny,D]=size(Y);
% D = size(Y, 2);
temp1 = exp(-sum((Y-V).^2,3)./(2*sigma2+sigma0));% Trick1:add small value to sigma to tolerant the turbulent/noise  
temp2 = (2*pi*sigma2)^(D/2)*(1-gamma)/(gamma*a);
P=temp1./(temp1+temp2);
P(~((I(:,:,1))))=0;P(~((I(:,:,2))))=0;
E=sum(sum(P.*sum((Y-V).^2,3)))/(2*sigma2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%
function V=smoothvfc(Y)
    V(:,:,1)=smoothn(Y(:,:,1),1000);
    V(:,:,2)=smoothn(Y(:,:,2),1000);
end

function demo()
    [x,y] = meshgrid(linspace(0,1,24));
    Vx = cos(2*pi*x+pi/2).*cos(2*pi*y);
    Vy = sin(2*pi*x+pi/2).*sin(2*pi*y);
    Vx = Vx + (0.01)*randn(24,24); % adding Gaussian noise
    Vy = Vy + (0.01)*randn(24,24); % adding Gaussian noise
    I = randperm(numel(Vx));
    Vx(I(1:60)) = (rand(60,1)-0.5)*5; % adding outliers
    Vy(I(1:60)) = (rand(60,1)-0.5)*5; % adding outliers
    Vx(I(61:100)) = NaN; % missing values
    Vy(I(61:100)) = NaN; % missing values
    VecFld = VFC(Vx,Vy,0.00);  %VFC Method
    Vx_VFC = VecFld.V(:,:,1);Vy_VFC = VecFld.V(:,:,2);VFC_Flag = VecFld.VFCIndex;% Output
    subplot(121), quiver(x,y,Vx,Vy,2.5), axis square
    title('Noisy velocity field')
    subplot(122), quiver(x,y,Vx_VFC,Vy_VFC), axis square
    title('Processed velocity field')
end
