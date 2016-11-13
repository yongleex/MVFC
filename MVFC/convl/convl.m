%% conventional method of post-processing of PIV data
% it consists of three steps:normalized median test + 3¡Á3 median
% replacement + 3¡Á3 smoothing kernel(with equal weights). The methods call
% the JPIV package 
% revised by Yong Lee @ 2015.09.02: add the outlier index output and
% smoothflag input
function [newu,newv,OutlierIndex] = convl(u0,v0,smoothflag)
if nargin<3
    smoothflag = true;% a default flag to do smooth
end
% default steps
dx=1;dy=1;
[Sizx,Sizy]=size(u0);% I will guarentee the size of u0 is equal to the size of v0
X = zeros(Sizx,Sizy);Y = zeros(Sizx,Sizy); OutlierIndex = zeros(Sizx,Sizy);
X = bsxfun(@plus,X,[1:dx:dx*Sizx]');
Y = bsxfun(@plus,Y,[1:dx:dx*Sizy]);
Input=[X(:),Y(:),u0(:),v0(:),ones(Sizx*Sizy,1)];
oh=jpiv2.PivData(Input);% object of java class handle
javaMethod('normalizedMedianTest',oh,0.1,2.0)%% normalized median test with 
javaMethod('replaceByMedian',oh,false,false)%% only replace the invalid vectors and exclude invalid vectors in the median calculation
if smoothflag
 javaMethod('smooth',oh,false)%% smooth with a 3¡Á3 averaging kernel
end


Output = javaMethod('getPivData',oh);%% get the results
newu = Output(:,3);
newv = Output(:,4);
OutlierIndex = Output(:,5);
newu = reshape(newu,[Sizx,Sizy]);
newv = reshape(newv,[Sizx,Sizy]);
OutlierIndex = reshape(OutlierIndex,[Sizx,Sizy]);
% OutlierIndex(OutlierIndex == -1) = 0; %change to the rule in this project (1: inlier; 0: outlier)
OutlierIndex = (OutlierIndex>=0);
end