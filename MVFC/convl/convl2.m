%% Conventional method of post-processing of PIV data
% it consists of three steps:normalized median test + 3¡Á3 median
% replacement + 3¡Á3 smoothing kernel(with equal weights). The author
% reprogrammed it with Matlab, the result is the same with it from the JPIV package.
% Yong Lee (leeyong@hust.edu.cn)
% 2015.10.18
% Input: the original velocity field u0, v0; the noise level (noiseLevel =
% 0.1 recommended); the threshold value (2, recommmended); window size (3
% or 5 recommend); Smooth or replace flags (set true if you want conduct the operation)

function [newu,newv,OutlierIndex] = convl2(u0,v0,noiseLevel, threshold,smoothflag,windowSize,ReplaceFlag)
[sx,sy] = size(u0);
newu = u0; newv = v0; OutlierIndex = (ones(sx,sy)>0);
%% Step 1.Median test
for j = 1:sy    % loop over every spacial point.
    for i = 1:sx
        for c = 1:2 % loop over velocity component u,v
            Neighbours = getNeighbours(); % Get the position of neighbours center at (i,j)
            sn = prod(size(Neighbours));  % The availabel neighbors number
            if sn>0
                if c == 1;Ux = u0;else;Ux = v0;end;
                if  OutlierIndex(i,j) == true
                    U  = Ux(Neighbours); U = U(:);
                    Temp = sort(U);
                    Um = (Temp(round((sn-1-0.25)/2)+1,:)+Temp(round((sn-1+0.25)/2)+1,:))/2; % Median value of the neighbours
                    
                    Res = U-Um;
                    Temp = sort(abs(Res));
                    Resm = (Temp(round((sn-1-0.25)/2)+1,:)+Temp(round((sn-1+0.25)/2)+1,:))/2;% Median value of the absolute value of residual
                    
                    OutlierIndex(i,j) = (abs(Ux(i,j)-Um)./(noiseLevel + Resm) <= threshold); % Normalized Median Test
                end
            else
                OutlierIndex(i,j)  = false;
            end
        end
    end
end

%% Step 2. Replacement of the outliers
if ReplaceFlag;replace(); end;

%% Step 3. Smooth operation
if smoothflag
    h    = fspecial('average',[3,3]);% Average filter only using the validated value
    Y    = filter2(h,double(OutlierIndex));
    newu = filter2(h,newu.*OutlierIndex)./Y; newu(isnan(newu)) = 0;
    newv = filter2(h,newv.*OutlierIndex)./Y; newv(isnan(newv)) = 0;
%     newu = medfilt2(newu,[3,3],'symmetric'); % Smooth the field with the median filter.
%     newv = medfilt2(newv,[3,3],'symmetric');
end

%% supporting function 1 get the neighbours
    function pos = getNeighbours()
        pos =[];
        ix1 = max(i-(windowSize-1)/2,1);ix2 = min(i+(windowSize-1)/2,sx);
        iy1 = max(j-(windowSize-1)/2,1);iy2 = min(j+(windowSize-1)/2,sy);
        for ix =ix1:ix2
            for jx= iy1:iy2
                if i == ix && j==jx;
                elseif OutlierIndex(ix,jx)
                    pos = [pos,(jx-1)*sx+ix];
                end
            end
        end
    end

%% supporting function 3 Replacement of the outlier with the median value
    function replace()
        windowSize = 3;
        for i = 1:sx
            for j = 1:sy     % loop over every spacial point.
                if OutlierIndex(i,j) == false
                    c =1;   pos = getNeighbours(); if prod(size(pos))>0; U = u0(pos); newu(i,j)  = median(U(:));else newu(i,j) = 0 ; end
                    c =2;   pos = getNeighbours(); if prod(size(pos))>0; U = v0(pos); newv(i,j)  = median(U(:));else newv(i,j) = 0 ; end
                end
            end
        end
    end

end