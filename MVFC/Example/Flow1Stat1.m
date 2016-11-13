% Case description
%-  Flow data  type: Synthetic Vortex cellular flow
%-  Noise dist type: Homogeneous Gaussian white noise with 0.01*Vmax standard deviation(STD)
%-  Outlier   Ratio: 10% clustered outlier + 0 missing data
%
%
%- Yong Lee 2016.08.10
%- leeyong@hust.edu.cn
%- This test program is based on the supplementary material of VFC, which
%can be download from my Website: <a
%   href="matlab:web('http://yong-lee.weebly.com/')">yong-lee.weebly.com/</a>

% Note that, this procedure is very time comsuming, about 10 hours for the
% 50 Mento-Calo simulations.
function Flow1Stat1()
%% Initialization
close all; clear; rng('default')
MonteNum = 2;% set 50 when run the fig in paper
% about 20 min for 2 MonteNum\

[L,H] = meshgrid(1:1:20,1:0.5:12);
sigma= 0.01; % 0.01: standard deviation of Gaussian noise
outlierRatio = 0.1;
Ns = numel(L);

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
        [x,y,Vx,Vy,u,v,OutlierIndex_Truth] = simulateVelocity(32,sigma,outlierRatio);
     
        %- Our Proposed Modified VFC (MVFC) Method
        VecFld = MVFC_test(Vx,Vy,0.001,L(j),H(j)); 
        Vx_MVFC = VecFld.V(:,:,1);Vy_MVFC = VecFld.V(:,:,2);    OutlierIndex_MVFC = VecFld.Index;% Output
        
        ssim_mvfc(i,j) = vssim(Vx_MVFC,Vy_MVFC,u,v);
        
        nrmse_mvfc(i,j) = nrmse(Vx_MVFC,Vy_MVFC,u,v);
        
        udc_mvfc(i,j) = L_udc(OutlierIndex_Truth,OutlierIndex_MVFC);
        
        odc_mvfc(i,j) = L_odc(OutlierIndex_Truth,OutlierIndex_MVFC);
    end
end
toc

%% Analysis the data and show in figs
% calculate the statistics variable
MEAN_NRMSE_MVFC = mean(nrmse_mvfc,1);    VAR_NRMSE_MVFC = var(nrmse_mvfc,1);          LB_NRMSE_MVFC=min(nrmse_mvfc,1);    UB_NRMSE_MVFC=max(nrmse_mvfc,1);

MEAN_SSIM_MVFC = mean(ssim_mvfc,1);      VAR_SSIM_MVFC = var(ssim_mvfc,1);            LB_SSIM_MVFC=min(ssim_mvfc,1);      UB_SSIM_MVFC=max(ssim_mvfc,1);

MEAN_ODC_MVFC = mean(odc_mvfc,1);        VAR_ODC_MVFC = var(odc_mvfc,1);              LB_ODC_MVFC = min(odc_mvfc,1);      UB_ODC_MVFC  = max(odc_mvfc,1);

MEAN_UDC_MVFC = mean(udc_mvfc,1);      VAR_UDC_MVFC = var(udc_mvfc,1);            LB_UDC_MVFC = min(udc_mvfc,1);      UB_UDC_MVFC = max(udc_mvfc,1);

%% Draw the results
% NRMSE
H1=figure();
mesh(L,H,reshape(MEAN_NRMSE_MVFC,size(L)));hold on; box on;
xlabel('\fontsize{14}L ');
ylabel('\fontsize{14}H ');
zlabel('\fontsize{14}NRMSE ');
set(gca,'fontsize',12); 

%SSIM
H2=figure();
mesh(L,H,reshape(MEAN_SSIM_MVFC,size(L)));hold on; box on;
xlabel('\fontsize{14}L ');
ylabel('\fontsize{14}H ');
zlabel('\fontsize{14}SSIM ');
set(gca,'fontsize',12); 

%UDC
H3=figure();
mesh(L,H,reshape(MEAN_UDC_MVFC,size(L)));hold on; box on;
xlabel('\fontsize{14}L ');
ylabel('\fontsize{14}H ');
zlabel('\fontsize{14}UDC ');
set(gca,'fontsize',12); 

%ODC
H4=figure();
mesh(L,H,reshape(MEAN_ODC_MVFC,size(L)));hold on; box on;
xlabel('\fontsize{14}L ');
ylabel('\fontsize{14}H ');
zlabel('\fontsize{14}ODC ');
set(gca,'fontsize',12); 

%TWDC
H5=figure();
mesh(L,H,reshape(MEAN_UDC_MVFC+MEAN_ODC_MVFC,size(L)));hold on; box on;
xlabel('\fontsize{14}L ');
ylabel('\fontsize{14}H ');
zlabel('\fontsize{14}UDC+ODC ');
set(gca,'fontsize',12); 

% %% Draw the results
% % NRMSE
% H1=figure();hold on;box on;xlim([0 60]);ylim([-0.01 1.01]);%title('\fontsize{12}NRMSE');
% plot(outlierRatio*100,MEAN_NRMSE_CON,'--d','LineWidth',2,'Color',[0.85,0.33,0.1], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.85,0.33,0.1]);
% % plot(outlierRatio*100,MEAN_NRMSE_CON1,'--h','LineWidth',2,'Color',[0.93,0.69,0.13], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.93,0.69,0.13]);
% plot(outlierRatio*100,MEAN_NRMSE_PPPIV,'-h','LineWidth',2,'Color',[0.93,0.69,0.13], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.93,0.69,0.13]);
% plot(outlierRatio*100,MEAN_NRMSE_VTM,'--p','LineWidth',2,'Color',[0.64,0.08,0.18], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.64,0.08,0.18]);
% plot(outlierRatio*100,MEAN_NRMSE_VFCS,'--s','LineWidth',2,'Color',[0.47,0.67,0.19] , 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.47,0.67,0.19]);
% plot(outlierRatio*100,MEAN_NRMSE_MVFC,'-o','LineWidth',2,'Color',[0.49,0.18,0.56], 'MarkerSize',8, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.49,0.18,0.56]);
% 
% H11=legend('\fontsize{14}NTM','\fontsize{14}DCT-PLS','\fontsize{14}VTM','\fontsize{14}VFC','\fontsize{14}MVFC',0);
% xlabel('\fontsize{14}Percentage of scattered outliers:% ');
% ylabel('\fontsize{14}NRMSE ');set(gca,'fontsize',12); set(H1,'position',[100 100 800 400]);
% 
% %SSIM
% H2=figure();hold on;box on;xlim([0 60]);ylim([-0.01 1.01]);%title('\fontsize{12}SSIM');
% plot(outlierRatio*100,MEAN_SSIM_CON,'--d','LineWidth',2,'Color',[0.85,0.33,0.1], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.85,0.33,0.1]);
% % plot(outlierRatio*100,MEAN_SSIM_CON1,'--h','LineWidth',2,'Color',[0.93,0.69,0.13], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.93,0.69,0.13]);
% plot(outlierRatio*100,MEAN_SSIM_PPPIV,'-h','LineWidth',2,'Color',[0.93,0.69,0.13], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.93,0.69,0.13]);
% plot(outlierRatio*100,MEAN_SSIM_VTM,'--p','LineWidth',2,'Color',[0.64,0.08,0.18], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.64,0.08,0.18]);
% plot(outlierRatio*100,MEAN_SSIM_VFCS,'--s','LineWidth',2,'Color',[0.47,0.67,0.19] , 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.47,0.67,0.19]);
% plot(outlierRatio*100,MEAN_SSIM_MVFC,'-o','LineWidth',2,'Color',[0.49,0.18,0.56], 'MarkerSize',8, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.49,0.18,0.56]);
% 
% H22=legend('\fontsize{14}CON','\fontsize{14}DCT-PLS','\fontsize{14}VTM','\fontsize{14}VFC','\fontsize{14}MVFC',0);
% xlabel('\fontsize{14}Percentage of scattered outliers:% ');
% ylabel('\fontsize{14}SSIM ');set(gca,'fontsize',12); set(H2,'position',[200 200 800 400]);
% 
% %UDC
% H3=figure();hold on;box on;xlim([0 60]);set(gca,'yscale','log');%ylim([-0.01 1.01]);%title('\fontsize{12}NRMSE');
% plot(outlierRatio*100,MEAN_UDC_CON+1,'--d','LineWidth',2,'Color',[0.85,0.33,0.1], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.85,0.33,0.1]);
% plot(outlierRatio*100,MEAN_UDC_VTM+1,'--p','LineWidth',2,'Color',[0.64,0.08,0.18], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.64,0.08,0.18]);
% plot(outlierRatio*100,MEAN_UDC_FADV+1,'--+','LineWidth',2,'Color',[0.3,0.75,0.93], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.3,0.75,0.93]);
% plot(outlierRatio*100,MEAN_UDC_VFC+1,'--s','LineWidth',2,'Color',[0.47,0.67,0.19] , 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.47,0.67,0.19]);
% plot(outlierRatio*100,MEAN_UDC_MVFC+1,'-o','LineWidth',2,'Color',[0.49,0.18,0.56], 'MarkerSize',8, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.49,0.18,0.56]);
% H11=legend('\fontsize{14}NMT','\fontsize{14}VTM','\fontsize{14}FADV','\fontsize{14}VFC','\fontsize{14}MVFC',0);
% xlabel('\fontsize{14}Percentage of scattered outliers:% ');
% ylabel('\fontsize{14}UDC+1 ');set(gca,'fontsize',12); set(H3,'position',[300 100 800 400]);
% 
% %ODC
% H4=figure();hold on;box on;xlim([0 60]);set(gca,'yscale','log');%ylim([-0.01 1.01]);%title('\fontsize{12}SSIM');
% plot(outlierRatio*100,MEAN_ODC_CON+1,'--d','LineWidth',2,'Color',[0.85,0.33,0.1], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.85,0.33,0.1]);
% plot(outlierRatio*100,MEAN_ODC_VTM+1,'--p','LineWidth',2,'Color',[0.64,0.08,0.18], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.64,0.08,0.18]);
% plot(outlierRatio*100,MEAN_ODC_FADV+1,'--+','LineWidth',2,'Color',[0.3,0.75,0.93], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.3,0.75,0.93]);
% plot(outlierRatio*100,MEAN_ODC_VFC+1,'--s','LineWidth',2,'Color',[0.47,0.67,0.19] , 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.47,0.67,0.19]);
% plot(outlierRatio*100,MEAN_ODC_MVFC+1,'-o','LineWidth',2,'Color',[0.49,0.18,0.56], 'MarkerSize',8, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.49,0.18,0.56]);
% H22=legend('\fontsize{14}NMT','\fontsize{14}VTM','\fontsize{14}FADV','\fontsize{14}VFC','\fontsize{14}MVFC',0);
% xlabel('\fontsize{14}Percentage of scattered outliers:% ');
% ylabel('\fontsize{14}ODC+1 ');set(gca,'fontsize',12); set(H4,'position',[400 200 800 400]);
% 
% %TWDC
% H5=figure();hold on;box on;xlim([0 60]);set(gca,'yscale','log');%ylim([-0.01 1.01]);%title('\fontsize{12}SSIM');
% plot(outlierRatio*100,MEAN_ODC_CON+MEAN_UDC_CON+1,'--d','LineWidth',2,'Color',[0.85,0.33,0.1], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.85,0.33,0.1]);
% plot(outlierRatio*100,MEAN_ODC_VTM+MEAN_UDC_VTM+1,'--p','LineWidth',2,'Color',[0.64,0.08,0.18], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.64,0.08,0.18]);
% plot(outlierRatio*100,MEAN_ODC_FADV+MEAN_UDC_FADV+1,'--+','LineWidth',2,'Color',[0.3,0.75,0.93], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.3,0.75,0.93]);
% plot(outlierRatio*100,MEAN_ODC_VFC+MEAN_UDC_VFC+1,'--s','LineWidth',2,'Color',[0.47,0.67,0.19] , 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.47,0.67,0.19]);
% plot(outlierRatio*100,MEAN_ODC_MVFC+MEAN_UDC_MVFC+1,'-o','LineWidth',2,'Color',[0.49,0.18,0.56], 'MarkerSize',8, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.49,0.18,0.56]);
% H22=legend('\fontsize{14}NMT','\fontsize{14}VTM','\fontsize{14}FADV','\fontsize{14}VFC','\fontsize{14}MVFC',0);
% xlabel('\fontsize{14}Percentage of scattered outliers:% ');
% ylabel('\fontsize{14}UDC+ODC+1');set(gca,'fontsize',12); set(H5,'position',[500 100 800 400]);

end

%% Supporting function
function [x,y,Vx,Vy,u,v,Idx]=simulateVelocity(N,sigma,outlierRatio)

U_max = 10;
Nc    = 6;
[x,y] =meshgrid(linspace(0,1,N));
Vx = U_max*cos(2*pi*x+pi/2).*cos(2*pi*y);
Vy = U_max*sin(2*pi*x+pi/2).*sin(2*pi*y);
u = Vx;v = Vy;
Idx = ones(size(u));
%   Corrupt the original flow
Vx = Vx + U_max*(sigma)*randn(size(Vx)); % adding Gaussian noise
Vy = Vy + U_max*(sigma)*randn(size(Vx)); % adding Gaussian noise

I = randperm(numel(Vx));
n = round(2*outlierRatio*numel(Vx)/Nc);% more seeding points
if n>0
    Index(1,:) = I(1:n);% ÿ���ĵ�һ��
    for i = 2:Nc
        Temp = randi([0,1],[1,n]);
        ix = randi([1,(i-1)]);
        Index(i,:) = Index(i-ix,:) +  Temp.*(2*randi([0,1],[1,n])-1)+(1-Temp).*(2*randi([0,1],[1,n])-1)*N;
        % ����Ƿ��������
        Haha = find(Index(i,:)<1 | Index(i,:)> N*N);
        Hehe = find(ismember(Index(i,:),Index(1:(i-1),:)) == 1);
        VarIndex = [Haha,Hehe];
        
        while length(VarIndex) > 0
            xx = length(VarIndex);
            Temp = randi([0,1],[1,xx]);
            ix = randi([1,(i-1)]);
            Index(i,VarIndex) = Index(i-ix,VarIndex) +  Temp.*(2*randi([0,1],[1,xx])-1)+(1-Temp).*(2*randi([0,1],[1,xx])-1)*N;
            % ����Ƿ��������
            Haha = find(Index(i,:)< 1 | Index(i,:)> N*N);
            Hehe = [];
            for k = VarIndex
                Hehe = [Hehe, find(ismember(Index(i,k),Index(1:(i-1),k)) == 1)];
            end
            VarIndex = [Haha,Hehe];
        end
    end
    
    Ix = 1; temp = 10;
    for i = 1:n
        if abs(numel(unique(Index(:,1:i)))./N^2-outlierRatio)<temp
            temp = abs(numel(unique(Index(:,1:Ix)))./N^2-outlierRatio);
            Ix = i;
        end
    end
    n =Ix;
    
    Vx(Index(1,1:n)) = (rand(1,n)-0.5)*2*U_max; % adding outliers
    Vy(Index(1,1:n)) = (rand(1,n)-0.5)*2*U_max; % adding outliers
    Idx(Index(1,1:n))= 0;
    for i = 1:n
        Vx(Index(2:Nc,i)) = Vx(Index(1,i))+(rand(Nc-1,1)-0.5)*0.2*Vx(Index(1,i)); % adding outliers
        Vy(Index(2:Nc,i)) = Vy(Index(1,i))+(rand(Nc-1,1)-0.5)*0.2*Vy(Index(1,i)); % adding outliers
        Idx(Index(2:Nc,i)) = 0;
    end
end
Idx = (Idx>0);

end