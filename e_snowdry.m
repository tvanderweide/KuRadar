% e_snowdry.m
% HPM 09/19/03
% this function calculates the real and complex part of the dielectric
%   constant of dry snow, given the density, frequency and temperature
%   After [Tiuri et al, 1984]
% INPUT: rho = dry density of snow [kg/m^3]
%          f = frequency [Hz] (default 10e9)
%          T = temperature [deg C] (default -10)
% OUTPUT: e_s = complex diel constant of the snow

function e_s = e_snowdry(rho,f,T)

if nargin == 1
    f=10e9; % [GHz] default frequency 
    T=-10; % [deg C] default temperature
elseif nargin == 2
    T=-10;  % [deg C] default temperature
end 
rho=rho/1000; % change units to [g/cc]
e_r=1+1.7*rho+0.7*rho.^2; % real part of diel const of dry snow
e_ice=1.59e6*(1./f+1.23e-14*sqrt(f)).*exp(0.036*T); % imaginary part of diel const of ice, from Tiuri, 1984 eq 6
e_i=(0.52*rho+0.62*rho.^2).*e_ice; % imaginary part of diel const of dry snow
e_s=e_r-e_i*j; % complex diel const of dry snow

