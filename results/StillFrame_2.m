%% Video
close all
clear all
clc

global epso sigma

%% Load %%
epso        = 8.85418781762e-12;

z_gnd       = load('z_gnd.dat');
NN          = load('N.dat');
dd          = load('d.dat');
SS          = load('Conductivity.dat');

N.r         = NN(1);
N.z         = NN(2);
N.t         = NN(3);

d.r         = dd(1);
d.z         = dd(2);
d.t         = dd(3);

sigma.za    = SS(1);
sigma.ra    = SS(2);
sigma.a     = .5e-14;
sigma.zb    = SS(4);
sigma.rb    = SS(5);
sigma.b     = 5e-14;
sigma.alpha = SS(7);

clear NN dd CC SS

% N.r         = 130;
% N.z         = 290;
% N.t         = 300001;
%
% d.r         = 500;
% d.z         = 250;
% d.t         = 1e-4;

t = (0:N.t-1)*d.t;
r = (0:N.r-1)*d.r;
z = (0:N.z-1)*d.z+z_gnd;


%% Expected Electric Field
E_th = load('Ec.dat');

SYM.r            = duplicate(r'*1e-3);                                     %_km
SYM.z            = z*1e-3;                                                 %_km

Extrema.E.min    = -E_th(1);
Extrema.E.max    = E_th(1);
Extrema.rhot.min = -2;
Extrema.rhot.max = 2;


dim = 2;


n=154951;

fprintf('n = %i\n',n);
Er   = load(['Er',num2str(dim),'d',num2str(n-1),'.dat']);
Ez   = load(['Ez',num2str(dim),'d',num2str(n-1),'.dat']);
rhos = 1e-9*load(['rhos',num2str(dim),'d',num2str(n-1),'.dat']);
rho  = 1e-9*load(['rho',num2str(dim),'d',num2str(n-1),'.dat']);

rhot             = rhos + rho;
Et               = Er.^2 + Ez.^2;
SYM.rhot         = duplicate(rhot);                            %_C/m^3
%     SYM.Et           = duplicate(1/2*log10(Et));                       %log(_V/m)
SYM.Et           = duplicate(Et.^.5);                              %_V/m

TMP(n) = max(max(rhot));
subplot(221);
imagesc(SYM.r,SYM.z,log10(abs(SYM.rhot'*1e9)));
set(gca,'XMinorTick','on','YMinorTick','on');
axis([-r(end)*1e-3 r(end)*1e-3 z(1)*1e-3 z(end)*1e-3])
xlabel('r (km)','fontsize',12)
ylabel('z (km)','fontsize',12)
axis xy
caxis([-6 1])
colorbar('location','WestOutside');
title('\rho_t (nC/m^3)','fontsize',12);

subplot(222);
imagesc(SYM.r,SYM.z,log10(abs(SYM.Et'*1e-5)));
set(gca,'XMinorTick','on','YMinorTick','on');
axis([-r(end)*1e-3 r(end)*1e-3 z(1)*1e-3 z(end)*1e-3])
xlabel('r (km)','fontsize',12)
ylabel('z (km)','fontsize',12)
axis xy
caxis([-6 1])
colorbar('location','EastOutside');
title('|E_t| (kV/cm)','fontsize',12);

subplot(223);
%     plot(log10(abs(rhos(1,:)*1e9)),z*1e-3, 'r-', log10(abs(rho(1,:)*1e9)),z*1e-3, 'b-', log10(abs(rhot(1,:)*1e9)),z*1e-3, 'g:');
plot(rhos(1,:)*1e9,z*1e-3, 'r-', rho(1,:)*1e9,z*1e-3, 'b-', rhot(1,:)*1e9,z*1e-3, 'g:');
%     axis([-6 1 z(1)*1e-3 z(end)*1e-3])
axis([Extrema.rhot.min Extrema.rhot.max z(1)*1e-3 z(end)*1e-3])
%     xlabel('log(|\rho|) (nC/m^3)','fontsize',12)
xlabel('\rho (nC/m^3)','fontsize',12)
ylabel('z (km)','fontsize',12)
set(gca,'XMinorTick','on','YMinorTick','on')
legend('\rho','\rho_s','\rho_t','location','NorthEast');
legend('boxoff')

subplot(224);
%     plot(log10(abs(Ez(1,:)*1e-5)),z*1e-3, 'r-', log10(abs(E_th*1e-5)),z*1e-3, 'g-.',log10(abs(32/2.16*E_th*1e-5)),z*1e-3, 'g-.');
%     set(gca,'XMinorTick','on','YMinorTick','on')
%     axis([-8 1 z(1)*1e-3 z(end)*1e-3])
%     xlabel('E_z (kV/cm)','fontsize',12)
%     ylabel('z (km)','fontsize',12)
%     text(1, z(end)*1e-3,['t = ', num2str(t(n),'%6.1f'), ' s'],'fontsize',12,'VerticalAlignment','top', 'HorizontalAlignment','right');

plot(Ez(1,:)*1e-5,z*1e-3, 'r-', E_th*1e-5,z*1e-3, 'g-.',-E_th*1e-5,z*1e-3, 'g-.');
set(gca,'XMinorTick','on','YMinorTick','on')
axis([-1.5 1.5 z(1)*1e-3 z(end)*1e-3])
xlabel('E_z (kV/cm)','fontsize',12)
ylabel('z (km)','fontsize',12)
text(1.5, z(end)*1e-3,['t = ', num2str(t(n),'%6.1f'), ' s'],'fontsize',12,'VerticalAlignment','top', 'HorizontalAlignment','right');

% hgexport(gcf,'BJ.eps');