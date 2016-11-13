function [newu,newv,y] = vtmedian(u,v,K)

% An implementation of the paper: variable threshold outlier identification
% in PIV data.
% Input: u,v the velocity data in gridded matrix
% Output: outlier position (inlier: y=1; outlier: y=0)
% Parameter K (optional)
% Yong Lee, 2015.09.01
% leeyong@hust.edu.cn

if nargin <3
    K = 1;
end
y = vtmedian1(u).*vtmedian1(v);
y = (y>0); % label the outliers

h = fspecial('gaussian', 5,4*sqrt(2));h=h/h(3,3);% alpha with equal weight

newu = filter2(h,u.*y,'same');
newv = filter2(h,v.*y,'same');
Temp = filter2(h,y,'same');
newu = newu./Temp;  newv = newv./Temp; % Replaced by the Gaussian Filter value without the outliers

newu(y>0) = u(y>0); newv(y>0) = v(y>0); % The original data

OutlierIndex = y;
h    = fspecial('average',[3,3]);% Average filter using the validated value
Y    = filter2(h,double(OutlierIndex));
newu = filter2(h,newu.*OutlierIndex)./Y;
newv = filter2(h,newv.*OutlierIndex)./Y;


    function Index = vtmedian1(X)
        % aggressive outlier rejection method: 3*3 block, the median 3 value are inliers
        [sx,sy]= size(X);
        Index = zeros(sx,sy); OutIndexTemp = zeros(sx,sy); ProField = X;
        NeighboursX1 = -1*ones(sx,sy);   NeighboursX1(1,:)  = 0;
        NeighboursX2 =    ones(sx,sy);   NeighboursX2(sx,:) = 0;
        NeighboursY1 =-1*ones(sx,sy) ;   NeighboursY1(:,1)  = 0;
        NeighboursY2 =    ones(sx,sy);   NeighboursY2(:,sy) = 0;
        
        for i = 1:sx
            for j=1:sy
                TempXBlock = X((i+NeighboursX1(i,j)):(i+NeighboursX2(i,j)),(j+NeighboursY1(i,j)):(j+NeighboursY2(i,j)));
                [B,I] = sort(TempXBlock(:));
                PositionFlag = abs(NeighboursX1(i,j))+abs(NeighboursX2(i,j))+abs(NeighboursY1(i,j))+abs(NeighboursY2(i,j));
                if PositionFlag == 2, i1 = 2;i2 =3; end
                if PositionFlag == 3, i1 = 3;i2 =4; end
                if PositionFlag == 4, i1 = 4;i2 =6; end
                
                if X(i,j)<B(i1) |X(i,j) >B(i2)
                    OutIndexTemp(i,j) = 0; % label as an outlier
                else
                    OutIndexTemp(i,j) = 1; % label as an inlier
                end
            end
        end
        
        % Generating the provisonal velocity field: 5*5 weighted Gaussian
        % Filter
        h = fspecial('gaussian', 5,4*sqrt(2));h=h/h(3,3);% alpha with equal weight
        Y1 = filter2(h,X.*OutIndexTemp,'same');
        Y2 = filter2(h,OutIndexTemp,'same');
        ProField = Y1./Y2;
        
        % Extract the 'expected' difference between a velocity and its
        % immediate neighbours;8 neighbours
        T1 = zeros(sx,sy);
        NeighboursX1 = -1*ones(sx,sy);   NeighboursX1(1,:)  = 0;
        NeighboursX2 =    ones(sx,sy);   NeighboursX2(sx,:) = 0;
        NeighboursY1 =-1*ones(sx,sy) ;   NeighboursY1(:,1)  = 0;
        NeighboursY2 =    ones(sx,sy);   NeighboursY2(:,sy) = 0;
        
        for i=1:sx
            for j=1:sy
                ProFieldblock = ProField((i+NeighboursX1(i,j)):(i+NeighboursX2(i,j)),(j+NeighboursY1(i,j)):(j+NeighboursY2(i,j)));
                N = prod(size(ProFieldblock))-1;
                T1(i,j) = sum(abs(ProFieldblock(:)-ProField(i,j)))/N+K;
            end
        end
        
        % This threshold is then filtered by a Gaussian kernel as follows:9*9
        % kernel
        h = fspecial('gaussian', 9,4*sqrt(2));% alpha with equal weight
        T = filter2(h,T1,'same')./filter2(h,ones(sx,sy),'same');
        
        % Using this variable threshold to do outlier validation in a
        % local-median method: 3*3 windows
        for i=1:sx
            for j=1:sy
                Xblock = X((i+NeighboursX1(i,j)):(i+NeighboursX2(i,j)),(j+NeighboursY1(i,j)):(j+NeighboursY2(i,j)));
                um = median(Xblock(:));
                if abs(X(i,j)-um)>T(i,j)
                    Index(i,j) = 0;% outliers
                else
                    Index(i,j) = 1;% inliers
                end
            end
        end
    end
end