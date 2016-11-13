function OutlierIndex = fadv(u,v,Sc)
% An implementation of the paper: Flow-adaptive data validation scheme in
% PIV (FADV).
% Input: u,v the velocity data in gridded matrix
% Output: Outlier Index (inlier: y=1; outlier: y=0)
% Parameter K (optional)
% Yong Lee, 2015.09.01
% leeyong@hust.edu.cn

%% check the input
if nargin <3
   Sc = 10;
end

%% Both components validation seperately
  OutlierIndex = fadv1(u).*fadv1(v);
  [~,Ru] = fadv1(u);[~,Rv] = fadv1(u);
 

%% supported sub function to do validation
    function [Index,RX] = fadv1(X)
   %%Step 1: an aggressive outlier Rejection£¨this step is quite similar to the Shinneeb,2004 work£©
    % aggressive outlier rejection method: 3*3 block, the median 3 value are inliers 
            RX = X;
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
        ProField = Y1./Y2; % this is the so-called provisional velocity field
        
    %% Step 2: Conservative algorithm
     % Get a square region centered at every suspect vector
        NeighboursX1 = -1*ones(sx,sy);   NeighboursX1(1,:)  = 0;
        NeighboursX2 =    ones(sx,sy);   NeighboursX2(sx,:) = 0;
        NeighboursY1 =-1*ones(sx,sy) ;   NeighboursY1(:,1)  = 0;
        NeighboursY2 =    ones(sx,sy);   NeighboursY2(:,sy) = 0;
        Index = OutIndexTemp;
        for i = 1:sx
            for j = 1:sy
                if OutIndexTemp(i,j) == 0  % for every suspect vector
                   Xblock = ProField((i+NeighboursX1(i,j)):(i+NeighboursX2(i,j)),(j+NeighboursY1(i,j)):(j+NeighboursY2(i,j)));%This is the square region 3*3
                   Xi = Xblock(1-NeighboursX1(i,j),:);
                   Xj = Xblock(:,1-NeighboursY1(i,j));
                   Xblock(1-NeighboursX1(i,j),:) = []; Xblock(:,1-NeighboursY1(i,j))=[]; Xaverage = mean(Xblock(:));% average value
                   XiR = L_linearReg(Xi); XjR = L_linearReg(Xj);% linear regression value
                   delta_rms_i = norm(XiR-Xaverage)/sqrt(prod(size(Xi)))+0.001;% A uniform flow should be exactly the same  if we do not add this small value 
                   delta_rms_j = norm(XjR-Xaverage)/sqrt(prod(size(Xj)))+0.001;
                   if (abs(X(i,j)-XiR(1-NeighboursY1(i,j))) < Sc*delta_rms_i) && (abs(X(i,j)-XjR(1-NeighboursX1(i,j))) < Sc*delta_rms_j) % inlier
                       Index(i,j) = 1;
                   else % outlier
                       Index(i,j) = 0; %
                       RX(i,j)    = ProField(i,j);
                   end
                   
                end
            end
        end
    end
%% A supported function to do linear regression
    function z = L_linearReg(w)
        [sx,sy] = size(w);
        if sx == 1
            w = w';
        end
        A = ones(max(size(w)),2); A(:,1)=1:max(size(w));
        z = A*(A'*A)^-1*A'*w;
        z = reshape(z,[sx,sy]);
    end
 
end
