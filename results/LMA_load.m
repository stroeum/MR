% LMA data loading
% global d Init Links N Nb z_gnd

dxyz              = load('d.dat')*1e-3;                                    %_km
Nxyz              = load('N.dat');                                         %_dimensionless
z_gnd             = load('z_gnd.dat')*1e-3;                                %_km
Layers.data       = load('ChargeLayers.dat')*1e-3;                         %_kC, _km

% Derive main parameters
Nb.Layers         = size(Layers.data);
Nb.Layers         = Nb.Layers(1);

N.r               = Nxyz(1);
N.z               = Nxyz(2);
N.t               = Nxyz(3);

d.r               = dxyz(1);                                               %_km
d.z               = dxyz(2);                                               %_km
d.t               = dxyz(3);                                               %_km

L.r               = (N.r-1)*d.r;                                           %_km
L.z               = (N.z-1)*d.z;                                           %_km
L.t               = (N.t-1)*d.t;                                           %_km
clear Nxyz dxyz InitPoint

r                    = (0:N.r-1)*d.r;
z                    = (0:N.z-1)*d.z+z_gnd;
t                    = (0:N.t-1)*d.t;

SS          = load('Conductivity.dat');
zc          = SS(1)*1e-3;
rc          = SS(2)*1e-3;
sigma.a     = 5e-14;
za          = SS(4)*1e-3;
alpha       = SS(7)*1e-3;

for ii=1:N.r
    for jj=1:N.z
        sigma.dat(ii,jj)   = sigma.a*exp((z(jj)+z_gnd)/za)*( 1 - (1-tanh((r(ii)-rc)/alpha))/2 * (1- tanh((z(jj)-z_gnd-zc)/alpha))/2);
    end
end

clear Nxyz dxyz SS