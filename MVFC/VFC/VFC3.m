function VecFld = VFC3(u,v,w)
% this is the non-optimized version to do 3D3C velocity field, details refer to VFC
% Yong Lee £¨LeeYong@hust.edu.cn£©
% 2015.07.08
%% input

gamma=0.90; lambda=1; theta=0.5; a=[]; MaxIter=20; ecr=10^-5; minP=10^-5;
% H=figure();quiver(X(:,:,1),X(:,:,2),Y(:,:,1),Y(:,:,2));

%% begin the code
fprintf('Starting Vector Field Correction Learning:\n');
[Nx,Ny,Nz]=size(u); D=3;

% Initialization
Vu=smoothvfc(u);Vv=smoothvfc(v);Vw=smoothvfc(w);
iter=1;  tecr=1; E=1;K=[1 1 1 ;1 -8 1;1 1 1]; 
I =  isfinite(u+v+w);%% missing value support
u(~I)=Vu(~I);v(~I)=Vv(~I);w(~I)=Vw(~I);% give the missing data a value


sigma2 = (sum((u(:)-Vu(:)).^2)+sum((v(:)-Vv(:)).^2)+sum((w(:)-Vw(:)).^2))/(sum(I(:)));a=5*(sigma2)^(D/2);
%%
newE = [];%
newE2 = [];%
while (iter<MaxIter) && (tecr > ecr) && (sigma2 > 1e-4) 
    % E-step. 
    E_old=E;
    [P, E]=get_P(u,v,w,Vu,Vv,Vw, sigma2 ,gamma, a,I);
    newE = [newE,E];%   

    fprintf('iterate: %dth, gamma: %f, sigma2=%f\n', iter, gamma,sigma2);

    % M-step. 
    % -Updata the function f 
    warning off all
    Vu = smoothn(u,P); Vv = smoothn(v,P);Vw = smoothn(w,P);
    % -Updata sigma2
    Sp=sum(P(:));
    sigma2=sum(sum(sum(P.*((u-Vu).^2+(v-Vv).^2+(w-Vw).^2))))/(Sp*D);
    % -Update gamma
    numcorr = length(find(P(:) > theta));
    gamma = numcorr/(sum(I(:)));
    if gamma > 0.95, gamma = 0.95; end
    if gamma < 0.05, gamma = 0.05; end
  
    iter=iter+1;
end
    warning on all
    Vu = smoothn(u,P,0.05); Vv = smoothn(v,P,0.05);Vw = smoothn(w,P,0.05);
%% Output
VecFld.u = Vu;
VecFld.v = Vv;
VecFld.w = Vw;
VecFld.P = P;
VecFld.VFCIndex = (P > theta);
disp('Removing outliers succesfully completed.');
end



%%%%%%%%%%%%%%%%%%%%%%%%
function [P, E]=get_P(u,v,w,Vu,Vv,Vw, sigma2 ,gamma, a,I)
% GET_P estimates the posterior probability and part of the energy.
u(~I)=0;v(~I)=0;w(~I)=0;Vu(~I)=0;Vv(~I)=0;Vw(~I)=0;
[Nx,Ny,Nz]=size(u);D=3;

temp1 = exp(-((u-Vu).^2+(v-Vv).^2+(w-Vw).^2)/(2*sigma2+0.1));
temp2 = (2*pi*sigma2)^(D/2)*(1-gamma)/(gamma*a);
P=temp1./(temp1+temp2);
P(~((I(:,:,1))))=0;P(~((I(:,:,2))))=0;
% E=sum(sum(P.*sum((Y-V).^2,3)))/(2*sigma2);%+sum(P(:))*log(sigma2)*D/2;
E=sum(sum(sum(P.*(((u-Vu).^2+(v-Vv).^2+(w-Vw).^2)))))/(2*sigma2);%+sum(P(:))*log(sigma2)*D/2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%
function V=smoothvfc(Y)
    V=smoothn(Y,10000);
end
