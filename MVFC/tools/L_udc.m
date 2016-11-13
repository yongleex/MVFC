function [Number,Ratio] = L_udc(Truth,Test)
%% Un-detected Number and Ratio
% If the vector is an outlier, but labeled as an inlier, we call it a
% un-detected vector,and let the Number plus 1. The Ratio means the total
% Number Ratio with respect to the total outlier number in Truth.
% Input
%     Truth: the ture indicator of the outliers,1 represent the inliers and
%     0 denote the outliers
%     Test: the test indicator
% Note: L means the code is designed by Yong Lee (Leeyong@hust.edu.cn)
%       2015.09.02
%
% An example:
% Truth = ones(10,1); Test = ones(10,1); Test([3,5,7])=0; Truth([3,6]) = 0;
% [UDP_Number,UDP_Ratio] = L_udp(Truth,Test)

N0 = sum(Truth(:) == 0);% the outlier total number
Temp = Test(find(Truth==0));% all the outliers have the labels in Temp
Number = sum(Temp(:) == 1);
Ratio  = Number/N0;
end