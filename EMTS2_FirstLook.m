% First look at Ku Radar Data
% Thomas Van Der Weide
% Radar by Adrian Tang et al, JPL and UCLA, as part of NASA IIP20

%% Initializing Variables
clear all; close all;
flow=15e9; % [Hz] start freq
fhigh=15.5e9; % [Hz] stop freq
BW=fhigh-flow; % [Hz] bandwidth
Tpl=67.55e-6; % [s] pulse length
Fs=122.88e6; % [Hz] sample rate
v=3.0e8; % [m/s] speed in air
N=2^15; % number of points in FFT
w=(0:N/2-1)/(N)*Fs; % frequencies sampled
d=0.5*w*Tpl/(BW)*v; % Distance in Air


%% Load a profile
D=load('P:\SnowDrones\Surveys\2024\2024-03-26_GrandMesa\Radar\Ku\data2.csv');
D=D'; % transpose due to collection in row vectors
[nr,nc]=size(D);

ix=1:nc; iy=1:nr; d2=d(1:nr);
figure(3);clf;
imagesc(ix,d2,D,[0 2e10]); colorbar
%imagesc(ix,d2,log(TCAL),[0 2e5]); colorbar
%imagesc(ix,d2,TCAL,[0 2e2]); colorbar

