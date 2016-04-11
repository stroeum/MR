close all
clear all
clc

Q = load('Q.dat');

%%
clc
% Q(1,:) = (0:length(Q)-1);%*.01;
set(gcf,'Units','normalized');
get(gcf,'OuterPosition');
set(gcf,'OuterPosition',[0 .5 1 .5]);

% set(gcf,'Units','normalized','OuterPosition',[0 0 1 1])
plot(Q(1,:),Q(2,:),'r--',Q(1,:),Q(4,:),'b',Q(1,:),Q(6,:),'r',Q(1,:),Q(7,:),'b--',Q(1,:),Q(2,:)+Q(4,:)+Q(6,:)+Q(7,:),'k-','LineWidth',1) 
legend('Q_{LP}','Q_N','Q_P','Q_{SC}','Q_\Sigma','Location','EastOutside')
legend('boxoff');
% axis([0 40 -100 70]);
xlabel('t (s)','fontsize',18)
ylabel('Q (C)','fontsize',18)
set(gca,'fontsize',18,'XMinorTick','on','YMinorTick','on')
% title('Evolution of the net charge in each layer','fontsize',18);

hgexport(gcf,'~/Desktop/MR/G_tmp_long_run/Q.eps');