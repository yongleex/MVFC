% Draw figure function

xlim([1,20]);ylim([1,12]);zlim([0,40])
view(125,35)
figure_FontSize = 20;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(get(gca,'ZLabel'),'FontSize',figure_FontSize);

set(findobj('FontSize',10),'FontSize',figure_FontSize);