%% Video
close all
clear all
beep  off
clc

% global epso sigma Layers

%% Load %%
epso        = 8.85418781762e-12;

LMA_load
set(gcf,'Units','normalized','OuterPosition',[.25 .5 .75 .5])
%% Expected Electric Field
E_th = load('Ec.dat');
Ez2D = load('Ez2d541000.dat');

Ez   = Ez2D(1,:);

SYM.r            = duplicate(r');                                          %_km
SYM.z            = z;                                                      %_km
SYM.sigma        = duplicate(sigma.dat);

Layers.Line.Style = '-';
Layers.Line.Width = 1;
Layers.Edge.Color = [[1 0 0];[0 0 1];[1 0 0];[0 0 1]];

hold 
imagesc(SYM.r,SYM.z,log10(SYM.sigma'))
L.x = 40;
L.z = 30;
% axis([-r(end) r(end) z(1) z(end)]);

set(gca,'XMinorTick','on','YMinorTick','on')
colormap(1-gray)
caxis([-14 -11])
YTick = get(colorbar,'YTick');
set(get(colorbar,'XLabel'),'String','\sigma (S/m)');
set(colorbar,'YTickLabel',10.^YTick);

LMA_charges

% rectangle('Position',[-L.x/2 z_gnd L.x L.z],'EdgeColor','w')

Extrema.E.min    = -E_th(1);
Extrema.E.max    = E_th(1);

plot(Ez*1e-5*10,z,'b',-E_th*1e-5*10,z,'g--',E_th*1e-5*10,z,'g--');

axis image
axis([-12.5 12.5 z_gnd z_gnd+21])
XTick = get(gca,'XTick');
set(gca,'XTickLabel',abs(XTick)');
xlabel('r (km), 10\times E_z (kV/cm)','fontsize',12)
ylabel('z (km)','fontsize',12)

box on
% legend('Ez','E_{th}');
% legend('Location','NorthOutside')
hgexport(gcf,'~/Desktop/BJ/SimGeometry_BJ.eps');