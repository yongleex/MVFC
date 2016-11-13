function test_VTM()
u = ones(10);
u(5,5) =0.5; v =u; % generate a synthetic velocity field;
u(1,1) = 1.2;
OulierIndex = vtmedian(u,v,0.1)% test the variable treshold local-median method with the K=0.1
end