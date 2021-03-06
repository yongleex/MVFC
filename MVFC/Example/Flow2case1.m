function Flow2case1()
%% Visualized Displayment on Experimental Fields (Turbulent flow): Fig 7 in the paper
% To assess the performance of the method using Turbulent flow data,
% :download from PIV Challenge
% Yong Lee 2015.5.25 leeyong@hust.edu.cn
% Revised by Yong Lee 2015.09.18 leeyong@hust.edu.cn

%Ref1(VFC) : Yong Lee (2015) A Robust Vector Field Correction Method via a Mixture Statistical Model of PIV signal.submited to Exp Fluids
%Ref2(DCT) : Garcia D (2011) A fast all-in-one method for automated post-processing of PIV data. Exp Fluids
%Ref3(CON) : Jerry Westerweel(2005).Universal outlier detection for PIV data
%Ref4(FADV): Liu, Z. L., Jia, L. F., Zheng, Y., & Zhang, Q. K. (2008). Flow-adaptive data validation scheme in PIV. Chem Eng Sci 63(1), 1-11.
%Ref5(VTM) : Shinneeb, A. M., Bugg, J. D., & Balachandar, R. (2004). Variable threshold outlier identification in PIV data. Meas Sci Technol, 15(9), 1722-1732. 

clear;close all;

%% load the Original Flow %  u,v,OutlierIndex_Truth
load turbulent_jet_withlabel; % Fig 7 in the paper
% load vortex_pairWithLabel;
%  load vortex_pair

% Forming x,y,Vx,Vy,OutlierIndex_Truth
[x,y] = meshgrid(1:size(u,2),1:size(u,1));x = (x-1)/(size(u,1)-1);y = (y-1)/(size(u,2)-1);
Vx = u; Vy = v; 

% Vx = Vx- min(Vx(:))+1; Vy = Vy- min(Vy(:))+1;
% Forming the ground truth field
[u0,v0] = GenerateRef(u,v,OutlierIndex_Truth);

Rho = sqrt(Vx.^2+Vy.^2);
Theta1 = asin(Vx./Rho)/pi;

%% Post-Processing by PPPIV and Ours
% [Vx_CON,Vy_CON,OutlierIndex_CON] = convl(Vx,Vy,1);% The Conventional
% method from JPIV with smoothing
% [Vx_CON1,Vy_CON1,OutlierIndex_CON1] = convl(Vx,Vy,0);% The Conventional
% method from JPIV without smoothing
tic
noiseLevel = 0.1; threshold = 2.0; smoothflag = true; windowSize =  5; ReplaceFlag = true;
[Vx_CON,Vy_CON,OutlierIndex_CON]=convl2(Vx,Vy,noiseLevel, threshold,smoothflag,windowSize,ReplaceFlag);% The Conventional method with smoothing
% [Vx_CON,Vy_CON,OutlierIndex_CON]=convl2(Rho,Theta1,noiseLevel, threshold,smoothflag,windowSize,ReplaceFlag);% The Conventional method with smoothing
toc
noiseLevel = 0.1; threshold =2.0; smoothflag = false; windowSize = 5; ReplaceFlag = true;
[Vx_CON1,Vy_CON1,OutlierIndex_CON1]=convl2(Vx,Vy,noiseLevel, threshold,smoothflag,windowSize,ReplaceFlag);% The Conventional method without smoothing
% [Vx_CON1,Vy_CON1,OutlierIndex_CON1]=convl2(Rho,Theta1,noiseLevel, threshold,smoothflag,windowSize,ReplaceFlag);% The Conventional method without smoothing

tic
[Vx_DCT,Vy_DCT] = pppiv(Vx,Vy);% The PPPIV (DCT-PLS) method
% [Vx_DCT,Vy_DCT] = pppiv(Rho,Theta1);% The PPPIV (DCT-PLS) method
toc

tic
[~,~,OutlierIndex_VTMedian] = vtmedian(Vx,Vy,1);% Variable threshold local median method (VTM)
% [~,~,OutlierIndex_VTMedian] = vtmedian(Rho,Theta1,1);% Variable threshold local median method (VTM)
toc

OutlierIndex_FADV     = fadv(Vx,Vy,25);   % Flow-adaptive Data Validation method (FADV)
% OutlierIndex_FADV     = fadv(Rho,Theta1,25);   % Flow-adaptive Data Validation method (FADV)

tic
VecFld = VFC(Vx,Vy,1);  %Our Method (VFC)
% VecFld = VFC(Rho,Theta1,1);  %Our Method (VFC)
Vx_VFC = VecFld.V(:,:,1);Vy_VFC = VecFld.V(:,:,2);OutlierIndex_VFC = VecFld.VFCIndex;% Output of VFC
toc

tic
VecFld = MVFC(Vx,Vy,0.001);  %Our MMVFC Method
toc
% VecFld = MVFC(Rho,Theta1,0.1);  %Our MMVFC Method
Vx_MVFC = VecFld.V(:,:,1);Vy_MVFC = VecFld.V(:,:,2);    OutlierIndex_MVFC = VecFld.Index;% Output
      
%% Assessment of the methods on Outlier Index
% Un-detected outlier count and over-detected outlier count
UO_OutierCount = [  L_udc(OutlierIndex_Truth,OutlierIndex_CON),L_odc(OutlierIndex_Truth,OutlierIndex_CON);
    L_udc(OutlierIndex_Truth,OutlierIndex_VTMedian),L_odc(OutlierIndex_Truth,OutlierIndex_VTMedian);
    L_udc(OutlierIndex_Truth,OutlierIndex_FADV),L_odc(OutlierIndex_Truth,OutlierIndex_FADV);
    L_udc(OutlierIndex_Truth,OutlierIndex_VFC),L_odc(OutlierIndex_Truth,OutlierIndex_VFC);
    L_udc(OutlierIndex_Truth,OutlierIndex_MVFC),L_odc(OutlierIndex_Truth,OutlierIndex_MVFC);];
fprintf('The Ground truth of Outlier Number:%d\n',sum(~OutlierIndex_Truth(:)))
fprintf('Overdetected Number:%d(NMT);%d(VTM);%d(FADV);%d(VFC);%d(MVFC)\n',L_odc(OutlierIndex_Truth,OutlierIndex_CON),...
    L_odc(OutlierIndex_Truth,OutlierIndex_VTMedian),L_odc(OutlierIndex_Truth,OutlierIndex_FADV),L_odc(OutlierIndex_Truth,OutlierIndex_VFC),L_odc(OutlierIndex_Truth,OutlierIndex_MVFC));
fprintf('Undetected   Number:%d(NMT);%d(VTM);%d(FADV);%d(VFC);%d(MVFC)\n',L_udc(OutlierIndex_Truth,OutlierIndex_CON),...
    L_udc(OutlierIndex_Truth,OutlierIndex_VTMedian),L_udc(OutlierIndex_Truth,OutlierIndex_FADV),L_udc(OutlierIndex_Truth,OutlierIndex_VFC),L_udc(OutlierIndex_Truth,OutlierIndex_MVFC));


%% Assessment of the methods on Restorated Velocity Field
%   Calculate  the SSIM and NRMSE value of each method.
s1 = vssim(Vx,Vy,u0,v0);            e1 = nrmse(Vx,Vy,u0,v0);            % Noisy field 
s2 = vssim(Vx_CON,Vy_CON,u0,v0);    e2 = nrmse(Vx_CON,Vy_CON,u0,v0);    % CON Method 
s3 = vssim(Vx_DCT,Vy_DCT,u0,v0);    e3 = nrmse(Vx_DCT,Vy_DCT,u0,v0);    % DCT-PLS Method
s4 = vssim(Vx_VFC,Vy_VFC,u0,v0);    e4 = nrmse(Vx_VFC,Vy_VFC,u0,v0);    % Our VFC Method
s5 = vssim(Vx_MVFC,Vy_MVFC,u0,v0);  e5 = nrmse(Vx_MVFC,Vy_MVFC,u0,v0);    % MVFC Method


str1=['Relative structure similarity:',num2str(s1),'. And the Relative NRMSE: ',num2str(e1)];
str2=['Relative structure similarity:',num2str(s2),'. And the Relative NRMSE: ',num2str(e2)];
str3=['Relative structure similarity:',num2str(s3),'. And the Relative NRMSE: ',num2str(e3)];
str4=['Relative structure similarity:',num2str(s4),'. And the Relative NRMSE: ',num2str(e4)];
str5=['Relative structure similarity:',num2str(s5),'. And the Relative NRMSE: ',num2str(e5)];


%% Display the results with flow fields displayment
% Original flow with outlier labeled in red
figure;scrsz = get(0,'ScreenSize');set(gcf,'Position',scrsz);
quiver(x(OutlierIndex_Truth),y(OutlierIndex_Truth),u(OutlierIndex_Truth),v(OutlierIndex_Truth));hold on;
quiver(x(~OutlierIndex_Truth),y(~OutlierIndex_Truth),u(~OutlierIndex_Truth),v(~OutlierIndex_Truth),'r');title('Original flow with outliers: Labeled by hand');

figure;scrsz = get(0,'ScreenSize');set(gcf,'Position',scrsz);
quiver(x,y,u,v);title('Original flow with outliers');
hold on; contour(x,y,sqrt(u.^2+v.^2),[0.4 1 2 3 4 5 6 7 8],'b-');
text(0,-0.1,str1);xlim([-0.05,1.05]);ylim([-0.25,1.05]);box on;set(gca,'ytick',[]);set(gca,'xtick',[]);

% Our "ground truth" flow
figure;scrsz = get(0,'ScreenSize');set(gcf,'Position',scrsz);
quiver(u0,v0);title('Original flow without outliers: Replaced by weighted average of 5*5 Gaussian Kernel and smoothed');

% Other Methods and the results
figure;scrsz = get(0,'ScreenSize');set(gcf,'Position',scrsz);
% subplot(2,2,1);quiver(x,y,u,v);title('Original flow with outliers');
% hold on; contour(x,y,sqrt(u.^2+v.^2),[0.4 1 2 3 4 5 6 7 8],'b-');
% text(0,-0.1,str1);xlim([-0.05,1.05]);ylim([-0.25,1.05]);box on;set(gca,'ytick',[]);set(gca,'xtick',[]);
subplot(2,2,1);quiver(x,y,Vx_CON,Vy_CON);title('NMT');
hold on; contour(x,y,sqrt(Vx_CON.^2+Vy_CON.^2),[0.4 1 2 3 4 5 6 7 8],'b-');
text(0,-0.1,str2);xlim([-0.05,1.05]);ylim([-0.25,1.05]);box on;set(gca,'ytick',[]);set(gca,'xtick',[]);
subplot(2,2,2);quiver(x,y,Vx_DCT,Vy_DCT);title('DCT-PLS');
hold on; contour(x,y,sqrt(Vx_DCT.^2+Vy_DCT.^2),[0.4 1 2 3 4 5 6 7 8],'b-');
text(0,-0.1,str3);xlim([-0.05,1.05]);ylim([-0.25,1.05]);box on;set(gca,'ytick',[]);set(gca,'xtick',[]);
subplot(2,2,3);quiver(x,y,Vx_VFC,Vy_VFC);title('VFC(Ours)');
hold on; contour(x,y,sqrt(Vx_VFC.^2+Vy_VFC.^2),[0.4 1 2 3 4 5 6 7 8],'b-');
text(0,-0.2,str4);
text(0,-0.1,'Note that: the contour plot in paper is done using Origin 9.1');xlim([-0.05,1.05]);ylim([-0.25,1.05]);box on;set(gca,'ytick',[]);set(gca,'xtick',[]);
subplot(2,2,4);quiver(x,y,Vx_MVFC,Vy_MVFC);title('MVFC');
hold on; contour(x,y,sqrt(Vx_MVFC.^2+Vy_MVFC.^2),[0.4 1 2 3 4 5 6 7 8],'b-');
text(0,-0.2,str5);
xlim([-0.05,1.05]);ylim([-0.25,1.05]);box on;set(gca,'ytick',[]);set(gca,'xtick',[]);

% The turbulent energy specturm
sx = max(size(log(turb_energy_spectrum(u,v)+1))); w =1:sx;
H = figure;grid off;hold on;box on;xlim([1-0.5 sx+0.5]);
xlabel('\fontsize{14} Wave Number'); ylabel('\fontsize{14} (Energy+1) / Wave Number');title('\fontsize{18}Turbulent energy spectrum comparison')
plot(w,(turb_energy_spectrum(u0,v0)+1),'-^','LineWidth',2,'Color',[0.0 0.0 1.0], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.0 0.0 1.0]);%original flow
plot(w,(turb_energy_spectrum(Vx,Vy)+1),'--v','LineWidth',2,'Color',[0,0.45,0.74], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0,0.45,0.74]);%corrupted flow
plot(w,(turb_energy_spectrum(Vx_CON1,Vy_CON1)+1),'--h','LineWidth',2,'Color',[0.93,0.69,0.13], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.93,0.69,0.13]);%conventional methods without smoothness
plot(w,(turb_energy_spectrum(Vx_CON,Vy_CON)+1),'--d','LineWidth',2,'Color',[0.85,0.33,0.1], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.85,0.33,0.1]);%conventional methods
plot(w,(turb_energy_spectrum(Vx_DCT,Vy_DCT)+1),'--s','LineWidth',2,'Color',[0.49,0.18,0.56], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.49,0.18,0.56]);%DCT-PLS flow
plot(w,(turb_energy_spectrum(Vx_VFC,Vy_VFC)+1),'-o','LineWidth',2,'Color',[0.47,0.67,0.19] , 'MarkerSize',8, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.47,0.67,0.19]);%Ours flow
plot(w,(turb_energy_spectrum(Vx_MVFC,Vy_MVFC)+1),'--h','LineWidth',2,'Color',[0.93,0.69,0.13], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.93,0.69,0.13]);%MVFC

H11 = legend('\fontsize{14}Ground truth','\fontsize{14}Flow with Outliers','\fontsize{14}NMT(Without smooth)','\fontsize{14}NMT(Smoothed)','\fontsize{14}DCT-PLS','\fontsize{14}VFC (Ours)','\fontsize{14}MVFC');
set(gca,'Yscale','log'); set(gca,'fontsize',12)
set(H,'position',[ 100 100 800 500]);

%% Show the results of un-detected and over-detected outlier count
L_drawBar(UO_OutierCount);
figure; spy(~OutlierIndex_Truth);title('truth')
figure; spy(~OutlierIndex_CON);title('NMT')
figure; spy(~OutlierIndex_VTMedian);title('VTM')
figure; spy(~OutlierIndex_VFC);title('VFC')
figure; spy(~OutlierIndex_MVFC);title('MVFC')


%% Output the results to Origin 
Vx_Noise = Vx - Vx_MVFC; Vx_Noise(~OutlierIndex_MVFC) = 0;
Vy_Noise = Vy - Vy_MVFC; Vy_Noise(~OutlierIndex_MVFC) = 0;

Fig10a = [x(:),y(:),x(:)+0.005*u(:),y(:)+0.005*v(:),OutlierIndex_Truth(:) ];
Fig10b = [x(:),y(:),x(:)+0.10*Vx_Noise(:),y(:)+0.10*Vy_Noise(:)];
Fig10c = [x(:),y(:),x(:)+0.003*Vx(:),y(:)+0.003*Vy(:),OutlierIndex_VFC(:)];
Fig10d = [x(:),y(:),x(:)+0.003*Vx_CON(:),y(:)+0.003*Vy_CON(:)];
Fig10e = [x(:),y(:),x(:)+0.003*Vx_DCT(:),y(:)+0.003*Vy_DCT(:)];
Fig10f = [x(:),y(:),x(:)+0.003*Vx_VFC(:),y(:)+0.003*Vy_VFC(:)];

Map1 = [x(OutlierIndex_Truth(:)<0.5),y(OutlierIndex_Truth(:)<0.5)];
Map2 = [x(OutlierIndex_VFC(:)<0.5),y(OutlierIndex_VFC(:)<0.5)];
Map3 = [x(OutlierIndex_MVFC(:)<0.5),y(OutlierIndex_MVFC(:)<0.5)];
Map4 = [x(OutlierIndex_CON(:)<0.5),y(OutlierIndex_CON(:)<0.5)];
Map5 = [x(OutlierIndex_VTMedian(:)<0.5),y(OutlierIndex_VTMedian(:)<0.5)];

Fig1c = [x(:),y(:),x(:)+0.005*u(:),y(:)+0.005*v(:)];
Fig7a = [x(:),y(:),x(:)+0.005*Vx(:),y(:)+0.005*Vy(:),OutlierIndex_Truth(:)];
xx = x(46:61,1:16);  yy = y(46:61,1:16); Vxx = Vx(46:61,1:16); Vyy = Vy(46:61,1:16);OutlierIndex_Truthy=OutlierIndex_Truth(46:61,1:16);
Fig7b = [xx(:),yy(:),xx(:)+0.05*Vxx(:),yy(:)+0.05*Vyy(:),OutlierIndex_Truthy(:)];

Fig7c = [x(:),y(:),x(:)+0.005*Vx_CON(:),y(:)+0.005*Vy_CON(:),sqrt(Vx_CON(:).^2+Vy_CON(:).^2)];
Fig7d = [x(:),y(:),x(:)+0.005*Vx_DCT(:),y(:)+0.005*Vy_DCT(:),sqrt(Vx_DCT(:).^2+Vy_DCT(:).^2)];
Fig7e = [x(:),y(:),x(:)+0.005*Vx_VFC(:),y(:)+0.005*Vy_VFC(:),sqrt(Vx_VFC(:).^2+Vy_VFC(:).^2)];
end


%% An interpolation function of Velocity Flield (u,v). If the Label is false(0), we replace it with the weighted average of the neighbours.
function [u0,v0] = GenerateRef(u,v,OutlierIndex_Truth)

h = fspecial('gaussian', 5,4*sqrt(2));% Gaussian Kernel
Y1 = filter2(h,u.*OutlierIndex_Truth,'same');
Y2 = filter2(h,OutlierIndex_Truth,'same');
Y = Y1./Y2;
u0 = u; u0(~OutlierIndex_Truth) = Y(~OutlierIndex_Truth);% outliers are replaced with the Gaussian kernel on component u

Y1 = filter2(h,v.*OutlierIndex_Truth,'same');
Y2 = filter2(h,OutlierIndex_Truth,'same');
Y = Y1./Y2;
v0 = v; v0(~OutlierIndex_Truth) = Y(~OutlierIndex_Truth);% outliers are replaced with the Gaussian kernel on component v

% u0 = medfilt2(u0,[3,3]);% filter the field with a median filter
% v0 = medfilt2(v0,[3,3]);
h    = fspecial('average',[3,3]);
u0 = filter2(h,u0.*OutlierIndex_Truth)./filter2(h,OutlierIndex_Truth);
v0 = filter2(h,v0.*OutlierIndex_Truth)./filter2(h,OutlierIndex_Truth);
end