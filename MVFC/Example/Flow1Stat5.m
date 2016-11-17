%% Case description
%-  Flow data  type: Synthetic Vortex cellular flow
%-  Noise dist type: Homogeneous Gaussian white noise with 0.01*Vmax standard deviation(STD)
%-  Outlier   Ratio: 10% Scattered outlier + 0 missing data
%-  The variable   : Vortex number varies from 2 to 7
%
%To assess the performance of the method respecting to vortex resolution  
% on Un/over-detected count. The reference methods are introducted in the file Readme.txt
%- Yong Lee 2016.08.10
%- leeyong@hust.edu.cn
%- This test program is based on the supplementary material of VFC, which
%can be download from my Website: <a
%   href="matlab:web('http://yong-lee.weebly.com/')">yong-lee.weebly.com/</a>

function Flow1Stat5()
%% Initialization
close all; clear; rng('default')
MonteNum =2;% set 100 when run the fig in paper
sigma= 0.01;  %0.1,0.01: standard deviation of Gaussian noise
outlierRatio = 0.1;
Nv = 2:7
Ns = max(size(Nv));

% Output:Structure similarity(SSIM) and normalized root of
% mean-squared-error(NRMSE)
odp_pppiv=[];odp_vfcs=[];odp_con=[];odp_vtm=[];odp_fadv=[];
udp_pppiv=[];udp_vfcs=[];udp_con=[];udp_vtm=[];odp_fadv=[];
%% Run the data
tic
for i=1:MonteNum
    for j=1:Ns
        %     parfor j=1:Ns
        %- Generate the flow field
        [x,y,Vx,Vy,u,v,OutlierIndex_Truth] = simulateVelocity(32,sigma,outlierRatio,Nv(j));
        
        %- Normalised Median Test
        noiseLevel = 0.1; threshold =2; smoothflag = true; windowSize = 5; ReplaceFlag = true;
        [Vx_CON,Vy_CON,OutlierIndex_CON]=convl2(Vx,Vy,noiseLevel, threshold,smoothflag,windowSize,ReplaceFlag);% The Conventional method proposed by Jerry Westerweel
        
        noiseLevel = 0.1; threshold =2; smoothflag = false; windowSize = 5; ReplaceFlag = true;
        [Vx_CON1,Vy_CON1,~]=convl2(Vx,Vy,noiseLevel, threshold,smoothflag,windowSize,ReplaceFlag);% The Conventional method proposed by Jerry Westerweel
        
        %- DCT-PLS method
        [Vx_DCT,Vy_DCT]=pppiv(Vx,Vy);%PPPIV method
        
        %- VTM method
        [~,~,OutlierIndex_VTMedian] = vtmedian(Vx,Vy,1);% Variable threshold local median method
        
        %- FADV method
        OutlierIndex_FADV     = fadv(Vx,Vy,10);   % Flow-adaptive Data Validation method
        
        %- VFC method 
        VecFld = VFC(Vx,Vy,1);  
        Vx_VFC = VecFld.V(:,:,1);Vy_VFC = VecFld.V(:,:,2);    OutlierIndex_VFC = VecFld.VFCIndex;% Output
               
        %- Our Proposed Modified VFC (MVFC) Method
        VecFld = MVFC(Vx,Vy,0.001);  
        Vx_MVFC = VecFld.V(:,:,1);Vy_MVFC = VecFld.V(:,:,2);    OutlierIndex_MVFC = VecFld.Index;% Output
        
        %% odp and udp
        udp_con(i,j) = L_udc(OutlierIndex_Truth,OutlierIndex_CON);
        odp_con(i,j) = L_odc(OutlierIndex_Truth,OutlierIndex_CON);
        udp_vtm(i,j) = L_udc(OutlierIndex_Truth,OutlierIndex_VTMedian);
        odp_vtm(i,j) = L_odc(OutlierIndex_Truth,OutlierIndex_VTMedian);
        udp_fadv(i,j) = L_udc(OutlierIndex_Truth,OutlierIndex_FADV);
        odp_fadv(i,j) = L_odc(OutlierIndex_Truth,OutlierIndex_FADV);
        udp_vfcs(i,j) = L_udc(OutlierIndex_Truth,OutlierIndex_VFC);
        odp_vfcs(i,j) = L_odc(OutlierIndex_Truth,OutlierIndex_VFC);
        udp_mvfc(i,j) = L_udc(OutlierIndex_Truth,OutlierIndex_MVFC);
        odp_mvfc(i,j) = L_odc(OutlierIndex_Truth,OutlierIndex_MVFC);
    end
end
udp_CON = mean(udp_con,1);udp_VTM = mean(udp_vtm,1);udp_FADV = mean(udp_fadv,1);udp_VFCS = mean(udp_vfcs,1);udp_MVFC = mean(udp_mvfc,1);
odp_CON = mean(odp_con,1);odp_VTM = mean(odp_vtm,1);odp_FADV = mean(odp_fadv,1);odp_VFCS = mean(odp_vfcs,1);odp_MVFC = mean(odp_mvfc,1);
for i = 1:Ns
    fprintf('Vortex cell:%d*%d\n',Nv(i));
    fprintf('Overdetected Number:\n%f(NMT);%f(VTM);%f(FADV);%f(VFC);%f(MVFC)\n',odp_CON(i),odp_VTM(i),odp_FADV(i),odp_VFCS(i),odp_MVFC(i));
    fprintf('Undetected   Number:\n%f(NMT);%f(VTM);%f(FADV);%f(VFC);%f(MVFC)\n',udp_CON(i),udp_VTM(i),udp_FADV(i),udp_VFCS(i),udp_MVFC(i));
end

% save flow1stat5Results


%% Assessment of the methods on Outlier Index
UO_OutierCount = [ udp_CON(1),odp_CON(1);udp_VTM(1),odp_VTM(1);udp_FADV(1),odp_FADV(1);udp_VFCS(1),odp_VFCS(1);udp_MVFC(1),odp_MVFC(1);0,0;
    udp_CON(2),odp_CON(2);udp_VTM(2),odp_VTM(2);udp_FADV(2),odp_FADV(2);udp_VFCS(2),odp_VFCS(2);udp_MVFC(2),odp_MVFC(2);0,0;
    udp_CON(3),odp_CON(3);udp_VTM(3),odp_VTM(3);udp_FADV(3),odp_FADV(3);udp_VFCS(3),odp_VFCS(3);udp_MVFC(3),odp_MVFC(3);0,0;
    udp_CON(4),odp_CON(4);udp_VTM(4),odp_VTM(4);udp_FADV(4),odp_FADV(4);udp_VFCS(4),odp_VFCS(4);udp_MVFC(4),odp_MVFC(4);0,0;
    udp_CON(5),odp_CON(5);udp_VTM(5),odp_VTM(5);udp_FADV(5),odp_FADV(5);udp_VFCS(5),odp_VFCS(5);udp_MVFC(5),odp_MVFC(5);0,0;
    udp_CON(6),odp_CON(6);udp_VTM(6),odp_VTM(6);udp_FADV(6),odp_FADV(6);udp_VFCS(6),odp_VFCS(6);udp_MVFC(6),odp_MVFC(6);
    ];
L_drawBar(UO_OutierCount);

end


%% Supporting function
function [x,y,Vx,Vy,u,v,OutlierIndex_Truth]=simulateVelocity(N,sigma,outlierRatio,Nv)
U_max =10;
[x,y] =meshgrid(linspace(0,1,N));
Vx = U_max*cos(Nv*pi*x+pi/2).*cos(Nv*pi*y);
Vy = U_max*sin(Nv*pi*x+pi/2).*sin(Nv*pi*y);
OutlierIndex_Truth = ones(size(Vx));
u = Vx;v = Vy;
%   Corrupt the original flow
Vx = Vx + U_max*(sigma)*randn(size(Vx)); % adding Gaussian noise
Vy = Vy + U_max*(sigma)*randn(size(Vx)); % adding Gaussian noise

I = randperm(numel(Vx));
n = round(outlierRatio*numel(Vx));% the outlier ratio =0.50
Vx(I(1:n)) = (rand(n,1)-0.5)*4*U_max; % adding outliers
Vy(I(1:n)) = (rand(n,1)-0.5)*4*U_max; % adding outliers
OutlierIndex_Truth(I(1:n)) = 0;
end
