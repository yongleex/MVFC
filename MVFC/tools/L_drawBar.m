function L_drawBar(Y)
%% this is a function to draw Count number of un/over-detected outliers with different methods
% Input: Y is a Nm*2 matrix, represents the results of the Nm kings of
% methods,and the first col represent the undetected outlier count number,and the seNMTd col represents the over-detected outlier number; 
Y(:,1) = Y(:,1)+ 0; 
H1 = figure; set(H1,'position',[100 100 800 400]);
X = 1:(prod(size(Y))/2);
bar(X,Y,1,'stacked');
% xlim([0.5,4.5]);ylim([1,1000]);title('Wrong detected Vector Number')

legend('UDN','ODN');set(gca,'XTick',1:1:35);
L = {'NMT','VTM','FADV','VFC','MVFC','','NMT','VTM','FADV','VFC','MVFC','','NMT','VTM','FADV','VFC','MVFC','','NMT','VTM','FADV','VFC','MVFC','','NMT','VTM','FADV','VFC','MVFC','','NMT','VTM','FADV','VFC','MVFC'};
set(gca,'XTickLabel',L);
% set(gca,'Yscale','log');
set(gca,'fontsize',12)
h = gca;
th=rotateticklabel(h, 45);%
% legend('UDC','ODC');
legend('UDC','ODC','location','NorthWest');
end