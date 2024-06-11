% Script for plotting the radar data
% Includes time and frequency domain plots and various ways of looking
% at the data
% Code from early on in the project
% Probably written by HP

%%
clear all;
close all;

% Load File
%skyCalibration=csvread('H:\RedPitaya_C\RadarData\22-05-04_FieldTest\sky_calibration.csv');
%skyCal= mean(skyCalibration);
D = load( 'P:\SnowDrones\Surveys\2024\2023-12-21_MCS\Radar\Data1.csv' );
% calc variables
Fs=122.88e6; % [Hz] sample frequency
Tpl=66e-6; % [s] pulse length
B=0.3e9; % [Hz] bandwidth
Ns=floor(Fs*Tpl);
Nt=floor(length(D)/Ns);
%D2=D-skyCal;
D2=D;

r.tdata=D2';
r.nfft=2^15;
r.Fs=Fs;
r.alpha=2.0;

r2=cal_psd_radar(r);
TWT=r2.w*Tpl/B;
d=TWT*3e8/2;

D3=highpass(r2.PDATA',5e5,Fs);
D4=envelope(D3);

% Plot "frequency" domain
figure(3);clf
imagesc(r2.PDATA); colorbar %r2.PDATA
% caxis([-70 0 ])
len = size(D2,1);
% axis([1 len 0 500])
set(gca,'LineWidth',2, 'Fontsize',14)
xlabel('trace #')
ylabel('distance in air [m]')



%% look at time domain
t=1/Fs:1/Fs:Ns/Fs;
t_cropped = t(1:len);
figure(2);clf;plot(t_cropped,D2(:,11:12))
set(gca,'LineWidth',2, 'Fontsize',14)
xlabel('time [s]')
ylabel('amplitude [V]')

%%
%Use this to get an idea of what the caxis range should be
figure(4);clf;plot(r2.PDATA(:,100:110))

%% frequency domain
% Plot "frequency" domain
figure(3);clf
imagesc(r2.PDATA); colorbar %r2.PDATA
% caxis([70 120])
len = size(D2,1);
axis([1 len 0 400])
set(gca,'LineWidth',2, 'Fontsize',14)
xlabel('trace #')
ylabel('distance in air [m]')
%% 

% INPUT: tdata = time domain matrix (from get_tdata)
%            nfft = number of points in FFT
%           Fs = sample frequency [Hz]
%        alpha = [5.3] KaiserBessel parameter
% OUTPUT: PDATA = power spectral density estimates [dB]
%           w = frequencies sampled [Hz]
% SNTX: obj = cal_psd_radar(obj)


