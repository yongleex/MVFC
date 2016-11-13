function test_FADV()
%% a test function to check the FADV code
%  Yong Lee (leeyong@hust.edu.cn)
%  2015.09.02
u = ones(10);
u(5,5) =0.9; v =u; % generate a synthetic velocity field;
u(1,1) =0.9;
OulierIndex = fadv(u,v,10)% test the Flow-adaptive data validation method with the Sc=10
end