%% Case description
%-  Flow data  type: Synthetic Vortex cellular flow
%-  Noise dist type: Homogeneous Gaussian white noise with 0.01*Vmax standard deviation(STD)
%-  Outlier   Ratio: 0-60% Scattered outlier +  0 cluster + 0 missing data
%
%To assess the performance of the method respecting to scatter outliers percentage effect. The reference methods are introducted in the file Readme.txt
%- Yong Lee 2016.08.10
%- leeyong@hust.edu.cn
%- This test program is based on the supplementary material of VFC, which
%can be download from my Website: <a
%   href="matlab:web('http://yong-lee.weebly.com/')">yong-lee.weebly.com/</a>

function Flow1Stat2()
%% Initialization
close all; clear; rng('default')
MonteNum = 2;% set 100 when run the fig in paper
sigma= 0.01; % 0.01: standard deviation of Gaussian noise
outlierRatio = [0:0.025:0.6];
Ns = max(size(outlierRatio));

% Output:Structure similarity(SSIM) and normalized root of
% mean-squared-error(NRMSE)
nrmse_pppiv=[];nrmse_vfcs=[];nrmse_con=[];
ssim_pppiv=[];ssim_vfcs=[];ssim_con=[];
odc_pppiv=[];odc_vfcs=[];odc_con=[];
udc_pppiv=[];udc_vfcs=[];udc_con=[];

%% Run the data
tic
for i=1:MonteNum
%     parfor j=1:Ns
    for j=1:Ns
        %- Generate the flow field
        [x,y,Vx,Vy,u,v,OutlierIndex_Truth] = simulateVelocity(32,sigma,outlierRatio(j));
        
        %- Normalised Median Test
        noiseLevel = 0.1; threshold =2; smoothflag = true; windowSize = 5; ReplaceFlag = true;
        [Vx_CON,Vy_CON,OutlierIndex_CON]=convl2(Vx,Vy,noiseLevel, threshold,smoothflag,windowSize,ReplaceFlag);% The Conventional method proposed by Jerry Westerweel
        
        %- DCT-PLS method
        [Vx_DCT,Vy_DCT]=pppiv(Vx,Vy);%PPPIV method
        
        %- VTM method
        [Vx_VTM,Vy_VTM,OutlierIndex_VTM] = vtmedian(Vx,Vy,1);% Variable threshold local median method
        
        %- FADV method
        OutlierIndex_FADV     = fadv(Vx,Vy,10);   % Flow-adaptive Data Validation method
        
        %- VFC method 
        VecFld = VFC(Vx,Vy,1.0); 
        Vx_VFC = VecFld.V(:,:,1);Vy_VFC = VecFld.V(:,:,2); OutlierIndex_VFC =VecFld.VFCIndex; % Output
        
        %- Our Proposed Modified VFC (MVFC) Method
        VecFld = MVFC(Vx,Vy,0.001); 
        Vx_MVFC = VecFld.V(:,:,1);Vy_MVFC = VecFld.V(:,:,2);    OutlierIndex_MVFC = VecFld.Index;% Output
        
        ssim_con(i,j) = vssim(Vx_CON,Vy_CON,u,v);        
        ssim_pppiv(i,j) = vssim(Vx_DCT,Vy_DCT,u,v);      ssim_vfcs(i,j) = vssim(Vx_VFC,Vy_VFC,u,v);
        ssim_vtm(i,j) = vssim(Vx_VTM,Vy_VTM,u,v);        ssim_mvfc(i,j) = vssim(Vx_MVFC,Vy_MVFC,u,v);
        
        nrmse_con(i,j) = nrmse(Vx_CON,Vy_CON,u,v);        
        nrmse_pppiv(i,j) = nrmse(Vx_DCT,Vy_DCT,u,v);      nrmse_vfcs(i,j) = nrmse(Vx_VFC,Vy_VFC,u,v);
        nrmse_vtm(i,j) = nrmse(Vx_VTM,Vy_VTM,u,v);        nrmse_mvfc(i,j) = nrmse(Vx_MVFC,Vy_MVFC,u,v);
        
        udc_con(i,j)  = L_udc(OutlierIndex_Truth,OutlierIndex_CON);        udc_vtm(i,j) = L_udc(OutlierIndex_Truth,OutlierIndex_VTM);
        udc_fadv(i,j) = L_udc(OutlierIndex_Truth,OutlierIndex_FADV);       udc_vfc(i,j) = L_udc(OutlierIndex_Truth,OutlierIndex_VFC);
        udc_mvfc(i,j) = L_udc(OutlierIndex_Truth,OutlierIndex_MVFC);
        
        odc_con(i,j)  = L_odc(OutlierIndex_Truth,OutlierIndex_CON);        odc_vtm(i,j) = L_odc(OutlierIndex_Truth,OutlierIndex_VTM);
        odc_fadv(i,j) = L_odc(OutlierIndex_Truth,OutlierIndex_FADV);       odc_vfc(i,j) = L_odc(OutlierIndex_Truth,OutlierIndex_VFC);
        odc_mvfc(i,j) = L_odc(OutlierIndex_Truth,OutlierIndex_MVFC);
    end
end
toc

%% Analysis the data and show in figs
% calculate the statistics variable
MEAN_NRMSE_CON = mean(nrmse_con,1);      VAR_NRMSE_CON = var(nrmse_con,1);            LB_NRMSE_CON=min(nrmse_con,1);      UB_NRMSE_CON=max(nrmse_con,1);
MEAN_NRMSE_PPPIV = mean(nrmse_pppiv,1);  VAR_NRMSE_PPPIV = var(nrmse_pppiv,1);        LB_NRMSE_PPPIV=min(nrmse_pppiv,1);  UB_NRMSE_PPPIV=max(nrmse_pppiv,1);
MEAN_NRMSE_VTM = mean(nrmse_vtm,1);      VAR_NRMSE_VTM  = var(nrmse_vtm,1);           LB_NRMSE_VTM  = min(nrmse_vtm,1);   UB_NRMSE_VTM =max(nrmse_vtm,1);
MEAN_NRMSE_VFCS = mean(nrmse_vfcs,1);    VAR_NRMSE_VFCS = var(nrmse_vfcs,1);          LB_NRMSE_VFCS=min(nrmse_vfcs,1);    UB_NRMSE_VFCS=max(nrmse_vfcs,1);
MEAN_NRMSE_MVFC = mean(nrmse_mvfc,1);    VAR_NRMSE_MVFC = var(nrmse_mvfc,1);          LB_NRMSE_MVFC=min(nrmse_mvfc,1);    UB_NRMSE_MVFC=max(nrmse_mvfc,1);

MEAN_SSIM_CON = mean(ssim_con,1);        VAR_SSIM_CON = var(ssim_con,1);              LB_SSIM_CON=min(ssim_con,1);        UB_SSIM_CON=max(ssim_con,1);
MEAN_SSIM_PPPIV = mean(ssim_pppiv,1);    VAR_SSIM_PPPIV = var(ssim_pppiv,1);          LB_SSIM_PPPIV=min(ssim_pppiv,1);    UB_SSIM_PPPIV=max(ssim_pppiv,1);
MEAN_SSIM_VTM  = mean(ssim_vtm,1);       VAR_SSIM_VTM = var(ssim_vtm,1);              LB_SSIM_VTM=min(ssim_vtm,1);        UB_SSIM_VTM=max(ssim_vtm,1);
MEAN_SSIM_VFCS = mean(ssim_vfcs,1);      VAR_SSIM_VFCS = var(ssim_vfcs,1);            LB_SSIM_VFCS=min(ssim_vfcs,1);      UB_SSIM_VFCS=max(ssim_vfcs,1);
MEAN_SSIM_MVFC = mean(ssim_mvfc,1);      VAR_SSIM_MVFC = var(ssim_mvfc,1);            LB_SSIM_MVFC=min(ssim_mvfc,1);      UB_SSIM_MVFC=max(ssim_mvfc,1);


MEAN_ODC_CON = mean(odc_con,1);          VAR_ODC_CON = var(odc_con,1);                LB_ODC_CON = min(odc_con,1);        UB_ODC_CON  = max(odc_con,1);
MEAN_ODC_VTM = mean(odc_vtm,1);          VAR_ODC_VTM = var(odc_vtm,1);                LB_ODC_VTM = min(odc_vtm,1);        UB_ODC_VTM  = max(odc_vtm,1);
MEAN_ODC_FADV= mean(odc_fadv,1);         VAR_ODC_FADV= var(odc_fadv,1);               LB_ODC_FADV= min(odc_fadv,1);       UB_ODC_FADV = max(odc_fadv,1);
MEAN_ODC_VFC = mean(odc_vfc,1);          VAR_ODC_VFC = var(odc_vfc,1);                LB_ODC_VFC = min(odc_vfc,1);        UB_ODC_VFC  = max(odc_vfc,1);
MEAN_ODC_MVFC = mean(odc_mvfc,1);        VAR_ODC_MVFC = var(odc_mvfc,1);              LB_ODC_MVFC = min(odc_mvfc,1);      UB_ODC_MVFC  = max(odc_mvfc,1);


MEAN_UDC_CON = mean(udc_con,1);        VAR_UDC_CON = var(udc_con,1);              LB_UDC_CON = min(udc_con,1);        UB_UDC_CON = max(udc_con,1);
MEAN_UDC_VTM = mean(udc_vtm,1);        VAR_UDC_VTM = var(udc_vtm,1);              LB_UDC_VTM = min(udc_vtm,1);        UB_UDC_VTM = max(udc_vtm,1);
MEAN_UDC_FADV= mean(udc_fadv,1);       VAR_UDC_FADV= var(udc_fadv,1);             LB_UDC_FADV= min(udc_fadv,1);       UB_UDC_FADV= max(udc_fadv,1);
MEAN_UDC_VFC = mean(udc_vfc,1);        VAR_UDC_VFC = var(udc_vfc,1);              LB_UDC_VFC = min(udc_vfc,1);        UB_UDC_VFC = max(udc_vfc,1);
MEAN_UDC_MVFC = mean(udc_mvfc,1);      VAR_UDC_MVFC = var(udc_mvfc,1);            LB_UDC_MVFC = min(udc_mvfc,1);      UB_UDC_MVFC = max(udc_mvfc,1);


%% Draw the results
% NRMSE
H1=figure();hold on;box on;xlim([0 60]);ylim([-0.01 1.01]);%title('\fontsize{12}NRMSE');
plot(outlierRatio*100,MEAN_NRMSE_CON,'--d','LineWidth',2,'Color',[0.85,0.33,0.1], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.85,0.33,0.1]);
% plot(outlierRatio*100,MEAN_NRMSE_CON1,'--h','LineWidth',2,'Color',[0.93,0.69,0.13], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.93,0.69,0.13]);
plot(outlierRatio*100,MEAN_NRMSE_PPPIV,'-h','LineWidth',2,'Color',[0.93,0.69,0.13], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.93,0.69,0.13]);
plot(outlierRatio*100,MEAN_NRMSE_VTM,'--p','LineWidth',2,'Color',[0.64,0.08,0.18], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.64,0.08,0.18]);
plot(outlierRatio*100,MEAN_NRMSE_VFCS,'--s','LineWidth',2,'Color',[0.47,0.67,0.19] , 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.47,0.67,0.19]);
plot(outlierRatio*100,MEAN_NRMSE_MVFC,'-o','LineWidth',2,'Color',[0.49,0.18,0.56], 'MarkerSize',8, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.49,0.18,0.56]);

H11=legend('\fontsize{14}NMT','\fontsize{14}DCT-PLS','\fontsize{14}VTM','\fontsize{14}VFC','\fontsize{14}MVFC',0);
xlabel('\fontsize{14}Percentage of scattered outliers:% ');
ylabel('\fontsize{14}NRMSE ');set(gca,'fontsize',12); set(H1,'position',[100 100 800 400]);

%SSIM
H2=figure();hold on;box on;xlim([0 60]);ylim([-0.01 1.01]);%title('\fontsize{12}SSIM');
plot(outlierRatio*100,MEAN_SSIM_CON,'--d','LineWidth',2,'Color',[0.85,0.33,0.1], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.85,0.33,0.1]);
% plot(outlierRatio*100,MEAN_SSIM_CON1,'--h','LineWidth',2,'Color',[0.93,0.69,0.13], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.93,0.69,0.13]);
plot(outlierRatio*100,MEAN_SSIM_PPPIV,'-h','LineWidth',2,'Color',[0.93,0.69,0.13], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.93,0.69,0.13]);
plot(outlierRatio*100,MEAN_SSIM_VTM,'--p','LineWidth',2,'Color',[0.64,0.08,0.18], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.64,0.08,0.18]);
plot(outlierRatio*100,MEAN_SSIM_VFCS,'--s','LineWidth',2,'Color',[0.47,0.67,0.19] , 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.47,0.67,0.19]);
plot(outlierRatio*100,MEAN_SSIM_MVFC,'-o','LineWidth',2,'Color',[0.49,0.18,0.56], 'MarkerSize',8, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.49,0.18,0.56]);

H22=legend('\fontsize{14}NMT','\fontsize{14}DCT-PLS','\fontsize{14}VTM','\fontsize{14}VFC','\fontsize{14}MVFC',0);
xlabel('\fontsize{14}Percentage of scattered outliers:% ');
ylabel('\fontsize{14}SSIM ');set(gca,'fontsize',12); set(H2,'position',[200 200 800 400]);

%UDC
H3=figure();hold on;box on;xlim([0 60]);set(gca,'yscale','log');%ylim([-0.01 1.01]);%title('\fontsize{12}NRMSE');
plot(outlierRatio*100,MEAN_UDC_CON+1,'--d','LineWidth',2,'Color',[0.85,0.33,0.1], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.85,0.33,0.1]);
plot(outlierRatio*100,MEAN_UDC_VTM+1,'--p','LineWidth',2,'Color',[0.64,0.08,0.18], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.64,0.08,0.18]);
plot(outlierRatio*100,MEAN_UDC_FADV+1,'--+','LineWidth',2,'Color',[0.3,0.75,0.93], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.3,0.75,0.93]);
plot(outlierRatio*100,MEAN_UDC_VFC+1,'--s','LineWidth',2,'Color',[0.47,0.67,0.19] , 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.47,0.67,0.19]);
plot(outlierRatio*100,MEAN_UDC_MVFC+1,'-o','LineWidth',2,'Color',[0.49,0.18,0.56], 'MarkerSize',8, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.49,0.18,0.56]);
H11=legend('\fontsize{14}NMT','\fontsize{14}VTM','\fontsize{14}FADV','\fontsize{14}VFC','\fontsize{14}MVFC',0);
xlabel('\fontsize{14}Percentage of scattered outliers:% ');
ylabel('\fontsize{14}UDC+1 ');set(gca,'fontsize',12); set(H3,'position',[300 100 800 400]);

%ODC
H4=figure();hold on;box on;xlim([0 60]);set(gca,'yscale','log');%ylim([-0.01 1.01]);%title('\fontsize{12}SSIM');
plot(outlierRatio*100,MEAN_ODC_CON+1,'--d','LineWidth',2,'Color',[0.85,0.33,0.1], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.85,0.33,0.1]);
plot(outlierRatio*100,MEAN_ODC_VTM+1,'--p','LineWidth',2,'Color',[0.64,0.08,0.18], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.64,0.08,0.18]);
plot(outlierRatio*100,MEAN_ODC_FADV+1,'--+','LineWidth',2,'Color',[0.3,0.75,0.93], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.3,0.75,0.93]);
plot(outlierRatio*100,MEAN_ODC_VFC+1,'--s','LineWidth',2,'Color',[0.47,0.67,0.19] , 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.47,0.67,0.19]);
plot(outlierRatio*100,MEAN_ODC_MVFC+1,'-o','LineWidth',2,'Color',[0.49,0.18,0.56], 'MarkerSize',8, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.49,0.18,0.56]);
H22=legend('\fontsize{14}NMT','\fontsize{14}VTM','\fontsize{14}FADV','\fontsize{14}VFC','\fontsize{14}MVFC',0);
xlabel('\fontsize{14}Percentage of scattered outliers:% ');
ylabel('\fontsize{14}ODC+1 ');set(gca,'fontsize',12); set(H4,'position',[400 200 800 400]);

%TWDC
H5=figure();hold on;box on;xlim([0 60]);set(gca,'yscale','log');%ylim([-0.01 1.01]);%title('\fontsize{12}SSIM');
plot(outlierRatio*100,MEAN_ODC_CON+MEAN_UDC_CON+1,'--d','LineWidth',2,'Color',[0.85,0.33,0.1], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.85,0.33,0.1]);
plot(outlierRatio*100,MEAN_ODC_VTM+MEAN_UDC_VTM+1,'--p','LineWidth',2,'Color',[0.64,0.08,0.18], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.64,0.08,0.18]);
plot(outlierRatio*100,MEAN_ODC_FADV+MEAN_UDC_FADV+1,'--+','LineWidth',2,'Color',[0.3,0.75,0.93], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.3,0.75,0.93]);
plot(outlierRatio*100,MEAN_ODC_VFC+MEAN_UDC_VFC+1,'--s','LineWidth',2,'Color',[0.47,0.67,0.19] , 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.47,0.67,0.19]);
plot(outlierRatio*100,MEAN_ODC_MVFC+MEAN_UDC_MVFC+1,'-o','LineWidth',2,'Color',[0.49,0.18,0.56], 'MarkerSize',8, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.49,0.18,0.56]);
H22=legend('\fontsize{14}NMT','\fontsize{14}VTM','\fontsize{14}FADV','\fontsize{14}VFC','\fontsize{14}MVFC',0);
xlabel('\fontsize{14}Percentage of scattered outliers:% ');
ylabel('\fontsize{14}UDC+ODC+1');set(gca,'fontsize',12); set(H5,'position',[500 100 800 400]);

end

%% Supporting function
function [x,y,Vx,Vy,u,v,Idx]=simulateVelocity(N,sigma,outlierRatio)
U_max = 10;
[x,y] =meshgrid(linspace(0,1,N));
Vx = U_max*cos(2*pi*x+pi/2).*cos(2*pi*y);
Vy = U_max*sin(2*pi*x+pi/2).*sin(2*pi*y);
u = Vx;v = Vy;
%   Corrupt the original flow
Vx = Vx + U_max*(sigma)*randn(size(Vx)); % adding Gaussian noise
Vy = Vy + U_max*(sigma)*randn(size(Vx)); % adding Gaussian noise

I = randperm(numel(Vx));
n = round(outlierRatio*numel(Vx));% the outlier ratio =0.50
Vx(I(1:n)) = (rand(n,1)-0.5)*4*U_max; % adding outliers
Vy(I(1:n)) = (rand(n,1)-0.5)*4*U_max; % adding outliers
Idx = ones(size(Vx)); Idx(I(1:n)) = 0; Idx = (Idx>0);% the ground truth label of outliers
Vx=Vx;
end