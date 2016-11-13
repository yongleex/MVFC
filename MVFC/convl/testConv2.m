function testConv2()
u0 = ones(5,8); v0 = ones(5,8) +5*rand(5,8); 
% rand('state',200);
%  u0(1:2,1:2) = rand(2);v0(1:2,1:2) = rand(2);
noiseLevel = 0.1; threshold =2.0; smoothflag = true; windowSize = 3; ReplaceFlag = true;

[Vx_CON,Vy_CON,OutlierIndex_CON1] = convl(u0,v0,0);% The Conventional method proposed by Jerry Westerweel
[newu,newv,OutlierIndex2] = convl2(u0,v0,noiseLevel, threshold,smoothflag,windowSize,ReplaceFlag);


OutlierIndex_CON1 %standard 1 
OutlierIndex2
OutlierIndex2-OutlierIndex_CON1

end