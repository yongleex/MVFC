%% Case description
%-  Flow data  type: Synthetic Vortex cellular flow
%-  Noise dist type: Non-Homogeneous Gaussian white noise with 0.06*Vmax standard deviation(STD)
%-  Outlier   Ratio: 15% Scattered outlier +  0 cluster + 0 missing data
%
% This is an example to illustrate the post-processing accuracy with the
%state of art methods. The reference methods are introducted in the file Readme.txt
%- Yong Lee 2016.08.10
%- leeyong@hust.edu.cn
%- This test program is based on the supplementary material of VFC, which
%can be download from my Website: <a
%   href="matlab:web('http://yong-lee.weebly.com/')">yong-lee.weebly.com/</a>

clear;close all;rng('default')
%% Generate the Original Flow
U_max = 10;
N = 2;
[x,y] = meshgrid(linspace(0,1,32));
Vx = U_max*cos(N*pi*x+pi/2).*cos(N*pi*y);
Vy = U_max*sin(N*pi*x+pi/2).*sin(N*pi*y);
u=Vx;v=Vy;

%%   Corrupt the original flow
sigma = 0.06;%standard deviation of the corrupted  Gaussiannoise, the value is 0.01*Vmax
OutRatio = 0.10;% Outlier ratio
rng(100);    Vx = Vx + U_max*(sigma)*randn(size(Vx)).*x; % adding Gaussian noise
rng(150);    Vy = Vy + U_max*(sigma)*randn(size(Vx)).*x; % adding Gaussian noise
Vx_Noise = Vx-u; Vy_Noise = Vy-v;

OutlierIndex_Truth = ones(size(Vx));
rng(200);    I = randperm(numel(Vx));
n = round(OutRatio*numel(Vx));% the outlier ratio =0.50
rng(200);    Vx(I(1:n)) =Vx(I(1:n))+ (1-0.0*x(I(1:n))).*(rand(1,n)-0.5)*4*U_max; % adding outliers
rng(300);    Vy(I(1:n)) =Vy(I(1:n))+ (1-0.0*x(I(1:n))).*(rand(1,n)-0.5)*4*U_max; % adding outliers
OutlierIndex_Truth(I(1:n)) = 0;

% % A cluster of outliers
% Vx(17:20,17:20) = U_max; Vy(17:20,17:20)=-0.5*U_max; OutlierIndex_Truth(17:20,17:20) = 0;
% 
% randn('state',600);    I = randperm(numel(Vx));
% n = round(2*numel(Vx)/10);
% Vx(I(1:n)) = NaN; % adding missing value
% Vy(I(1:n)) = NaN; % adding missing value


%% Post-Processing by 5 methods
%- Normalised Median Test
noiseLevel = 0.1; threshold =2; smoothflag = true; windowSize = 5; ReplaceFlag = true;
[Vx_CON,Vy_CON,OutlierIndex_CON]=convl2(Vx,Vy,noiseLevel, threshold,smoothflag,windowSize,ReplaceFlag);% The Conventional method proposed by Jerry Westerweel

noiseLevel = 0.1; threshold =2; smoothflag = false; windowSize = 5; ReplaceFlag = true;
[Vx_CON1,Vy_CON1,OutlierIndex_CON1]=convl2(Vx,Vy,noiseLevel, threshold,smoothflag,windowSize,ReplaceFlag);% The Conventional method proposed by Jerry Westerweel

%- DCT-PLS method
[Vx_DCT,Vy_DCT]=pppiv(Vx,Vy);%PPPIV method

%- VTM method
[Vx_VTM,Vy_VTM,OutlierIndex_VTMedian] = vtmedian(Vx,Vy,1);% Variable threshold local median method

%- FADV method
OutlierIndex_FADV     = fadv(Vx,Vy,10);   % Flow-adaptive Data Validation method

%- VFC method 
VecFld = VFC(Vx,Vy,1);  %Our VFC Method
Vx_VFC = VecFld.V(:,:,1);Vy_VFC = VecFld.V(:,:,2);    OutlierIndex_VFC = VecFld.VFCIndex;% Output

%- Our Proposed Modified VFC (MVFC) Method
VecFld = MVFC(Vx,Vy,0.001);  %Our MMVFC Method
Vx_MVFC = VecFld.V(:,:,1);Vy_MVFC = VecFld.V(:,:,2);    OutlierIndex_MVFC = VecFld.Index;% Output

%% Assessment of the methods on Outlier Index
UO_OutierCount = [  L_udc(OutlierIndex_Truth,OutlierIndex_CON),L_odc(OutlierIndex_Truth,OutlierIndex_CON);
    L_udc(OutlierIndex_Truth,OutlierIndex_VTMedian),L_odc(OutlierIndex_Truth,OutlierIndex_VTMedian);
    L_udc(OutlierIndex_Truth,OutlierIndex_FADV),L_odc(OutlierIndex_Truth,OutlierIndex_FADV);
    L_udc(OutlierIndex_Truth,OutlierIndex_VFC),L_odc(OutlierIndex_Truth,OutlierIndex_VFC);
    L_udc(OutlierIndex_Truth,OutlierIndex_MVFC),L_odc(OutlierIndex_Truth,OutlierIndex_MVFC);
    ];

fprintf('Overdetected Number:%d(NMT/CON);%d(VTM);%d(FADV);%d(VFC);%d(MVFC)\n',L_odc(OutlierIndex_Truth,OutlierIndex_CON),...
    L_odc(OutlierIndex_Truth,OutlierIndex_VTMedian),L_odc(OutlierIndex_Truth,OutlierIndex_FADV),L_odc(OutlierIndex_Truth,OutlierIndex_VFC),L_odc(OutlierIndex_Truth,OutlierIndex_MVFC));
fprintf('Undetected   Number:%d(NMT/CON);%d(VTM);%d(FADV);%d(VFC);%d(MVFC)\n',L_udc(OutlierIndex_Truth,OutlierIndex_CON),...
    L_udc(OutlierIndex_Truth,OutlierIndex_VTMedian),L_udc(OutlierIndex_Truth,OutlierIndex_FADV),L_udc(OutlierIndex_Truth,OutlierIndex_VFC),L_udc(OutlierIndex_Truth,OutlierIndex_MVFC));


%% Assessment of the methods on Restored Velocity fields
%   calculate the index
s0 = vssim(Vx_CON,Vy_CON,u,v);  % CON Method
s1 = vssim(Vx,Vy,u,v);          % Noisy field
s2 = vssim(Vx_DCT,Vy_DCT,u,v);  % DCT-PLS Method
s3 = vssim(Vx_VFC,Vy_VFC,u,v);  % Our VFC Method
s4 = vssim(Vx_MVFC,Vy_MVFC,u,v);  % Our VFC Method


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


%% Display the results
% Original flow with outlier labeled in red
figure;quiver(x,y,u,v);title('Original flow without degradation');
text(0,-0.1,str);xlim([-0.05,1.05]);ylim([-0.25,1.05]);box on;set(gca,'ytick',[]);set(gca,'xtick',[]);

% Corrupted flow field
figure;quiver(x,y,Vx,Vy);title('Original flow corrupted with noise and outliers');
text(0,-0.2,str1);text(0,-0.1,['Noist STD:',num2str(sigma),'.Outlier Ratio:',num2str(OutRatio)]);
xlim([-0.05,1.05]);ylim([-0.25,1.05]);box on;set(gca,'ytick',[]);set(gca,'xtick',[]);

% The corresponding noise distribution
figure;quiver(x,y,Vx_Noise,Vy_Noise);title('The corresponding noise distribution');
xlim([-0.05,1.05]);ylim([-0.25,1.05]);box on;set(gca,'ytick',[]);set(gca,'xtick',[]);


% corrupted flow and processed flow
figure;scrsz = get(0,'ScreenSize');set(gcf,'Position',scrsz);
subplot(2,2,1);quiver(x,y,Vx_CON,Vy_CON);title('NMT/CON');
text(0,-0.1,str0);xlim([-0.05,1.05]);ylim([-0.25,1.05]);box on;set(gca,'ytick',[]);set(gca,'xtick',[]);

subplot(2,2,2);quiver(x,y,Vx_DCT,Vy_DCT);title('PPPIV/(DCT-PLS)');
text(0,-0.1,str2);xlim([-0.05,1.05]);ylim([-0.25,1.05]);box on;set(gca,'ytick',[]);set(gca,'xtick',[]);

subplot(2,2,3);quiver(x,y,Vx_VFC,Vy_VFC);title('VFC');
text(0,-0.1,str3);xlim([-0.05,1.05]);ylim([-0.25,1.05]);box on;set(gca,'ytick',[]);set(gca,'xtick',[]);

subplot(2,2,4);quiver(x,y,Vx_MVFC,Vy_MVFC);title('MVFC');
text(0,-0.2,str4);xlim([-0.05,1.05]);ylim([-0.25,1.05]);box on;set(gca,'ytick',[]);set(gca,'xtick',[]);

% The turbulent energy specturm
sx = max(size(log(turb_energy_spectrum(u,v)+1)));
w =1:sx;
H = figure;grid off;hold on;box on;xlim([1-0.5 sx+0.5]); xlabel('\fontsize{14} Wave Number'); ylabel('\fontsize{14} (Energy+1) / Wave Number')
title('\fontsize{18}Turbulent energy spectrum comparison')
plot(w,(turb_energy_spectrum(u,v)+1),'-^','LineWidth',2,'Color',[0.0 0.0 1.0], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.0 0.0 1.0]);%original flow
plot(w,(turb_energy_spectrum(Vx,Vy)+1),'--v','LineWidth',2,'Color',[0,0.45,0.74], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0,0.45,0.74]);%corrupted flow
plot(w,(turb_energy_spectrum(Vx_CON,Vy_CON)+1),'--d','LineWidth',2,'Color',[0.85,0.33,0.1], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.85,0.33,0.1]);%conventional methods
plot(w,(turb_energy_spectrum(Vx_DCT,Vy_DCT)+1),'--s','LineWidth',2,'Color',[0.49,0.18,0.56], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.49,0.18,0.56]);%DCT-PLS flow
plot(w,(turb_energy_spectrum(Vx_VFC,Vy_VFC)+1),'-o','LineWidth',2,'Color',[0.47,0.67,0.19] , 'MarkerSize',8, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.47,0.67,0.19]);%VFC
plot(w,(turb_energy_spectrum(Vx_MVFC,Vy_MVFC)+1),'--h','LineWidth',2,'Color',[0.93,0.69,0.13], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.93,0.69,0.13]);%MVFC

H11 = legend('\fontsize{14}Original Flow','\fontsize{14}Flow with Outliers','\fontsize{14}NMT/CON(Smoothed)','\fontsize{14}DCT-PLS','\fontsize{14}VFC','\fontsize{14}MVFC');
set(gca,'Yscale','log'); set(gca,'fontsize',12)
set(H,'position',[ 100 100 800 500]);

%% Show the results of un-detected and over-detected outlier count
L_drawBar(UO_OutierCount);
figure; spy(~flipud(OutlierIndex_Truth));title('truth')
figure; spy(~flipud(OutlierIndex_CON));title('NMT')
figure; spy(~flipud(OutlierIndex_VTMedian));title('VTM')
figure; spy(~flipud(OutlierIndex_VFC));title('VFC')
figure; spy(~flipud(OutlierIndex_MVFC));title('MVFC')

%% output for origin plot
Fig1a = [x(:),y(:),x(:)+0.003*u(:),y(:)+0.003*v(:)];
Fig2a = [x(:),y(:),x(:)+0.003*Vx(:),y(:)+0.003*Vy(:),OutlierIndex_Truth(:)];
Fig2b = [x(:),y(:),x(:)+0.02*Vx_CON(:),y(:)+0.02*Vy_CON(:)];
Fig2c = [x(:),y(:),x(:)+0.02*Vx_DCT(:),y(:)+0.02*Vy_DCT(:)];
Fig2d = [x(:),y(:),x(:)+0.02*Vx_VFC(:),y(:)+0.02*Vy_VFC(:)];

Map1 = [x(OutlierIndex_Truth(:)<0.5),y(OutlierIndex_Truth(:)<0.5)];
Map2 = [x(OutlierIndex_VFC(:)<0.5),y(OutlierIndex_VFC(:)<0.5)];
Map3 = [x(OutlierIndex_MVFC(:)<0.5),y(OutlierIndex_MVFC(:)<0.5)];


%     ÿ�ַ������ı�ʾ
% Original Flow ʵ������� ,'-^','LineWidth',2,'Color',[0.0 0.0 1.0], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.0 0.0 1.0]);%ori flow
% ��������      ��������� ,'--v','LineWidth',2,'Color',[0,0.45,0.74], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0,0.45,0.74]);%corrupted flow
% CON           ��������   ,'--d','LineWidth',2,'Color',[0.85,0.33,0.1], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.85,0.33,0.1]);%con1
%               ��������� ,'--h','LineWidth',2,'Color',[0.93,0.69,0.13], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.93,0.69,0.13]);%con2
% DCT-PLS       �������� ,'--s','LineWidth',2,'Color',[0.49,0.18,0.56], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.49,0.18,0.56]);%DCT-PLS flow
% VFC           ʵ��Բ��   ,'-o','LineWidth',2,'Color',[0.47,0.67,0.19] , 'MarkerSize',8, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.47,0.67,0.19] );%Ours flow
% FADV          ����+��    ,'--+','LineWidth',2,'Color',[0.3,0.75,0.93], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.3,0.75,0.93]);%FADV
% VTM           ��������� ,'--p','LineWidth',2,'Color',[0.64,0.08,0.18], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.64,0.08,0.18]);%VTM
