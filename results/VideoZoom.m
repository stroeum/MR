%% Video
close all
clear all
clc

global epso sigma

%% Load %%
epso        = 8.85418781762e-12;

NN          = load('N.dat');
dd          = load('d.dat');
SS          = load('Conductivity.dat');
z_gnd       = load('z_gnd.dat');

N.r         = NN(1);
N.z         = NN(2);
N.t         = NN(3);

d.r         = dd(1);
d.z         = dd(2);
d.t         = dd(3);

sigma.za    = SS(1);
sigma.ra    = SS(2);
sigma.a     = 0*1.5e-14;
sigma.zb    = SS(4);
sigma.rb    = SS(5);
sigma.b     = 5e-14;
sigma.alpha = SS(7);

clear NN dd SS

t = (0:N.t-1)*d.t;
r = (0:N.r-1)*d.r;
z = (0:N.z-1)*d.z+z_gnd;

%% Expected Electric Field
E_th = load('Ec.dat');

SYM.r            = duplicate(r'*1e-3);                                     %_km
SYM.z            = z*1e-3;                                                 %_km

imageNB          = 1;                                                      % Image counter
set(gcf,'Units','normalized','OuterPosition',[0.5 .25 .5 .75]);
set(gcf,'Color',[1 1 1]);
Movie(10)        = getframe(gcf);

dim = 2;
for n = 0:2000:N.t-1
    clf
    %         fprintf('n = %i\n',n);
    Er   = load(['Er',num2str(dim),'d',num2str(n),'.dat']);
    Ez   = load(['Ez',num2str(dim),'d',num2str(n),'.dat']);
    rhos = 1e-9*load(['rhos',num2str(dim),'d',num2str(n),'.dat']);
    rho  = 1e-9*load(['rho',num2str(dim),'d',num2str(n),'.dat']);

    rhot             = rhos + rho;
    Et               = (Er.^2 + Ez.^2).^(1/2);
    % SYM.Et           = duplicate(1/2*log10(Et*(1e-5)^2));              %log(_kV/cm)
    SYM.Et           = duplicate(Et*1e-5);              %log(_kV/cm)
    SYM.rhot         = duplicate(rhot*1e9);                            %_nC/m^3

    subplot(2,2,1:2);
    imagescSgnLog(SYM.r,SYM.z,SYM.rhot,-5,1);
    set(gca,'FontSize',18)
    set(gca,'XMinorTick','on','YMinorTick','on')
    title('\rho_t=\rho_s+\rho_f (nC/m^3)','fontsize',18);
    axis([-r(end)*1e-3 r(end)*1e-3 z(1)*1e-3 z(end)*1e-3])
    xlabel('r (km)','fontsize',18)
    ylabel('z (km)','fontsize',18)    
    XTick = get(gca,'XTick');
    set(gca,'XTickLabel',abs(XTick)');

    subplot(2,2,3);
    plot(rhos(1,:)*1e9,z*1e-3, 'r-', rho(1,:)*1e9,z*1e-3, 'b-', rhot(1,:)*1e9,z*1e-3, 'g--');
    set(gca,'XMinorTick','on','YMinorTick','on')
    axis([-2.5 2.5 z(1)*1e-3 (z(1)+20e3)*1e-3])
    xlabel('\rho (nC/m^3)','fontsize',18)
    ylabel('z (km)','fontsize',18)
    set(gca,'FontSize',18)

    legend('\rho_s','\rho_f','\rho_t','location','SouthEast');
    legend('boxoff')

    %         subplot(224);
    %         plot(log10(abs(Ez(1,:))),z*1e-3, 'r-', log10(abs(E_th)),z*1e-3, 'b--');
    %         set(gca,'XMinorTick','on','YMinorTick','on')
    %         axis([Extrema.E.min Extrema.E.max z(1)*1e-3 z(end)*1e-3])
    %         xlabel('log(|E_z|) (kV/cm)','fontsize',10)
    %         ylabel('z (km)','fontsize',10)
    %         text(Extrema.E.max, z(end)*1e-3,['t = ', num2str(t(n+1),'%6.1f'), ' s'],'fontsize',10,'VerticalAlignment','top', 'HorizontalAlignment','right');

    subplot(2,2,4);
    plot(Ez(1,:)*1e-5,z*1e-3, 'r-', E_th*1e-5,z*1e-3, 'g--',-E_th*1e-5,z*1e-3, 'g--');
    set(gca,'XMinorTick','on','YMinorTick','on')
    axis([-max(E_th)*1e-5 max(E_th)*1e-5 z(1)*1e-3 (z(1)+20e3)*1e-3])
    xlabel('E_z (kV/cm)','fontsize',18)
    ylabel('z (km)','fontsize',18)
    set(gca,'FontSize',18)
    text(max(E_th)*1e-5,(z(1)+20e3)*1e-3,['t = ', num2str(t(n+1),'%6.1f'), ' s'],'fontsize',18,'VerticalAlignment','top', 'HorizontalAlignment','right');

    Movie(imageNB)   = getframe(gcf);
    imageNB          = imageNB+1;

end

%% Animate
% Record = input('Record the movie? (1: yes, else: no)\n>> ');
Record = 1;
if (Record == 1)
    movie2avi(Movie,'Source.avi','quality',100);
else
    return
end
