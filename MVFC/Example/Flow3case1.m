%% Intuitive Visual Comparison
%This is an example to illustrate the post-processing accuracy with
%comparision with PPPIV(Garcia,2011).
% Yong Lee 2015.5.25
% leeyong@hust.edu.cn
%Ref1: Yong Lee (2015) A Robust Vector Field Correction Method via a Mixture Statistical Model of PIV signal.submited to Exp Fluids
%Ref2: Garcia D (2011) A fast all-in-one method for automated post-processing of PIV data. Exp Fluids
%Ref3: Jerry Westerweel(2005).Universal outlier detection for PIV data

clear;close all;
%% Generate the Original Flow : A source flow
U_max =10;r0 = 0.4;
[x,y] = meshgrid(linspace(-1,1,32));
R2=(x.^2+y.^2);
Source = (R2<r0^2);
R2(R2<r0^2) = r0^2;
Vx = U_max*r0*x./R2/2;
Vy = U_max*r0*y./R2/2;
u=Vx;v=Vy;

%%   Corrupt the original flow
sigma=0.01;%standard deviation of the corrupted Gaussian noise, the value is 0.01*Vmax
OutRatio =0.50;
randn('state',100);    Vx = Vx + U_max*(sigma)*randn(size(Vx)); % adding Gaussian noise
randn('state',150);    Vy = Vy + U_max*(sigma)*randn(size(Vx)); % adding Gaussian noise

rand('state',200);    I = randperm(numel(Vx));
n = round(OutRatio*numel(Vx));% the outlier ratio =0.50
rand('state',200);    Vx(I(1:n)) = (rand(n,1)-0.5)*2*U_max; % adding outliers
rand('state',300);    Vy(I(1:n)) = (rand(n,1)-0.5)*2*U_max; % adding outliers

OutlierIndex_Truth = ones(size(Vx));% Outlier labels
OutlierIndex_Truth(I(1:n)) = 0;

%     randn('state',600);    I = randperm(numel(Vx));
%     n = round(2*numel(Vx)/10);
%     Vx(I(1:n)) = NaN; % adding missing value
%     Vy(I(1:n)) = NaN; % adding missing value


%% Post-Processing by PPPIV and Ours
%
%     [Vx_CON,Vy_CON,OutlierIndex_CON]=convl(Vx,Vy,1);% The Conventional method xproposed by Jerry Westerweel
%     [Vx_CON1,Vy_CON1]=convl(Vx,Vy,0);% The Conventional method xproposed by Jerry Westerweel
noiseLevel = 0.1; threshold =2; smoothflag = true; windowSize = 5; ReplaceFlag = true;
[Vx_CON,Vy_CON,OutlierIndex_CON]=convl2(Vx,Vy,noiseLevel, threshold,smoothflag,windowSize,ReplaceFlag);% The Conventional method proposed by Jerry Westerweel

noiseLevel = 0.1; threshold =2; smoothflag = false; windowSize = 5; ReplaceFlag = true;
[Vx_CON1,Vy_CON1,~]=convl2(Vx,Vy,noiseLevel, threshold,smoothflag,windowSize,ReplaceFlag);% The Conventional method proposed by Jerry Westerweel

[Vx_DCT,Vy_DCT]=pppiv(Vx,Vy);%PPPIV method
[~,~,OutlierIndex_VTMedian] = vtmedian(Vx,Vy,1);% Variable threshold local median method
OutlierIndex_FADV     = fadv(Vx,Vy,10);   % Flow-adaptive Data Validation method
%     prod(size(Vx))-sum(OutIndex_VTMedian(:))
VecFld = VFC(Vx,Vy,1);  %Our Method
Vx_VFC = VecFld.V(:,:,1);Vy_VFC = VecFld.V(:,:,2);    OutlierIndex_VFC = VecFld.VFCIndex;% Output

VecFld = MVFC(Vx,Vy,0.001);  %MVFC
Vx_MVFC = VecFld.V(:,:,1);Vy_MVFC = VecFld.V(:,:,2);    OutlierIndex_MVFC = VecFld.Index;% Output

%% Assessment of the methods on Outlier Index
UO_OutierCount = [  L_udc(OutlierIndex_Truth,OutlierIndex_CON),L_odc(OutlierIndex_Truth,OutlierIndex_CON);
    L_udc(OutlierIndex_Truth,OutlierIndex_VTMedian),L_odc(OutlierIndex_Truth,OutlierIndex_VTMedian);
    L_udc(OutlierIndex_Truth,OutlierIndex_FADV),L_odc(OutlierIndex_Truth,OutlierIndex_FADV);
    L_udc(OutlierIndex_Truth,OutlierIndex_VFC),L_odc(OutlierIndex_Truth,OutlierIndex_VFC);
    L_udc(OutlierIndex_Truth,OutlierIndex_MVFC),L_odc(OutlierIndex_Truth,OutlierIndex_MVFC);];
L_drawBar(UO_OutierCount);

fprintf('Overdetected Number:%d(NMT/CON);%d(VTM);%d(FADV);%d(VFC);%d(VMFC)\n',L_odc(OutlierIndex_Truth,OutlierIndex_CON),...
    L_odc(OutlierIndex_Truth,OutlierIndex_VTMedian),L_odc(OutlierIndex_Truth,OutlierIndex_FADV),L_odc(OutlierIndex_Truth,OutlierIndex_VFC),L_odc(OutlierIndex_Truth,OutlierIndex_MVFC));
fprintf('Undetected   Number:%d(NMT/CON);%d(VTM);%d(FADV);%d(VFC);%d(MVFC)\n',L_udc(OutlierIndex_Truth,OutlierIndex_CON),...
    L_udc(OutlierIndex_Truth,OutlierIndex_VTMedian),L_udc(OutlierIndex_Truth,OutlierIndex_FADV),L_udc(OutlierIndex_Truth,OutlierIndex_VFC),L_odc(OutlierIndex_Truth,OutlierIndex_MVFC));


%% Assessment of the methods
%   calculate the index
s0 = vssim(Vx_CON,Vy_CON,u,v);
s1 = vssim(Vx,Vy,u,v);
s2 = vssim(Vx_DCT,Vy_DCT,u,v);
s3 = vssim(Vx_VFC,Vy_VFC,u,v);
s4 = vssim(Vx_MVFC,Vy_MVFC,u,v);


e0 = nrmse(Vx_CON,Vy_CON,u,v);
e1 = nrmse(Vx,Vy,u,v);
e2 = nrmse(Vx_DCT,Vy_DCT,u,v);
e3 = nrmse(Vx_VFC,Vy_VFC,u,v);
e4 = nrmse(Vx_MVFC,Vy_MVFC,u,v);

str=['structure similarity:',num2str(1),'. And the NRMSE: ',num2str(0)];
str0=['structure similarity:',num2str(s0),'. And the NRMSE: ',num2str(e0)];
str1=['structure similarity:',num2str(s1),'. And the NRMSE: ',num2str(e1)];
str2=['structure similarity:',num2str(s2),'. And the NRMSE: ',num2str(e2)];
str3=['structure similarity:',num2str(s3),'. And the NRMSE: ',num2str(e3)];
str4=['structure similarity:',num2str(s4),'. And the NRMSE: ',num2str(e4)];


%% show the results
% original flow
figure;quiver(x,y,u,v);title('Original flow without degradation');
text(-1,-1.1,str);xlim([-1.05,1.05]);ylim([-1.25,1.05]);box on;set(gca,'ytick',[]);set(gca,'xtick',[]);
text(-1,-1.2,'A uniform source in the circle');
rectangle('Position',[0-r0 0-r0 2*r0 2*r0],'Curvature',[1,1]);axis equal;
% corrupted flow and processed flow
figure;quiver(x,y,Vx,Vy);title('Original flow corrupted with noise and outliers');axis equal
text(-1,-1.2,str1);text(0,-1.1,['Noist STD:',num2str(sigma),'.Outlier Ratio:',num2str(OutRatio)]);
xlim([-1.05,1.05]);ylim([-1.25,1.05]);box on;set(gca,'ytick',[]);set(gca,'xtick',[]);

figure;scrsz = get(0,'ScreenSize');set(gcf,'Position',scrsz);
subplot(2,2,1);quiver(x,y,Vx_CON,Vy_CON);title('NMT/CON');axis equal
text(-1,-1.1,str0);xlim([-1.05,1.05]);ylim([-1.25,1.05]);box on;set(gca,'ytick',[]);set(gca,'xtick',[]);
% subplot(2,2,1);quiver(x,y,Vx,Vy);title('Original flow corrupted with noise and outliers');axis equal
% text(-1,-1.2,str1);text(0,-1.1,['Noist STD:',num2str(sigma),'.Outlier Ratio:',num2str(OutRatio)]);
% xlim([-1.05,1.05]);ylim([-1.25,1.05]);box on;set(gca,'ytick',[]);set(gca,'xtick',[]);
subplot(2,2,2);quiver(x,y,Vx_DCT,Vy_DCT);title('DCT-PLS');axis equal
text(-1,-1.1,str2);xlim([-1.05,1.05]);ylim([-1.25,1.05]);box on;set(gca,'ytick',[]);set(gca,'xtick',[]);
subplot(2,2,3);quiver(x,y,Vx_VFC,Vy_VFC);title('VFC');axis equal
text(-1,-1.1,str3);xlim([-1.05,1.05]);ylim([-1.25,1.05]);box on;set(gca,'ytick',[]);set(gca,'xtick',[]);
subplot(2,2,4);quiver(x,y,Vx_MVFC,Vy_MVFC);title('MVFC');axis equal
text(-1,-1.1,str3);xlim([-1.05,1.05]);ylim([-1.25,1.05]);box on;set(gca,'ytick',[]);set(gca,'xtick',[]);

% the turbulent energy specturm
sx = max(size(log(turb_energy_spectrum(u,v)+1)));
w =1:sx;
H = figure;grid off;hold on;box on;xlim([1-0.5 sx+0.5]); xlabel('\fontsize{14} Wave Number'); ylabel('\fontsize{14} (Energy+1) / Wave Number')
%     title('\fontsize{18}Turbulent energy spectrum comparison')
plot(w,(turb_energy_spectrum(u,v)+1),'-^','LineWidth',2,'Color',[0.0 0.0 1.0], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.0 0.0 1.0]);%original flow
plot(w,(turb_energy_spectrum(Vx,Vy)+1),'--v','LineWidth',2,'Color',[0.75 0.75 0], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.75 0.75 0]);%corrupted flow
plot(w,(turb_energy_spectrum(Vx_CON,Vy_CON)+1),'--d','LineWidth',2,'Color',[0.75 0 0.75], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.75 0 0.75]);%conventional methods
plot(w,(turb_energy_spectrum(Vx_DCT,Vy_DCT)+1),'--s','LineWidth',2,'Color',[0.85 0.33 0.1], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.85 0.33 0.1]);%DCT-PLS flow
plot(w,(turb_energy_spectrum(Vx_VFC,Vy_VFC)+1),'-o','LineWidth',2,'Color',[0.0 0.5 0], 'MarkerSize',8, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.0 0.5 0]);%Ours flow
plot(w,(turb_energy_spectrum(Vx_MVFC,Vy_MVFC)+1),'-h','LineWidth',2,'Color',[0.93,0.69,0.13], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.93,0.69,0.13]);%conventional methods without smoothness

H11 = legend('\fontsize{14}Original Flow','\fontsize{14}Flow with Outliers','\fontsize{14}NMT/CON','\fontsize{14}DCT-PLS','\fontsize{14}VFC','\fontsize{14}MVFC',0);

set(gca,'Yscale','log'); set(gca,'fontsize',12)
set(H,'position',[ 100 100 800 500]);

%% output for origin plot
Fig1b = [x(:),y(:),x(:)+0.2*u(:),y(:)+0.2*v(:),Source(:)];
Fig4a = [x(:),y(:),x(:)+0.2*Vx(:),y(:)+0.2*Vy(:),OutlierIndex_VFC(:)];
Fig4b = [x(:),y(:),x(:)+0.2*Vx_CON(:),y(:)+0.2*Vy_CON(:)];
Fig4c = [x(:),y(:),x(:)+0.2*Vx_DCT(:),y(:)+0.2*Vy_DCT(:)];
Fig4d = [x(:),y(:),x(:)+0.2*Vx_VFC(:),y(:)+0.2*Vy_VFC(:)];
