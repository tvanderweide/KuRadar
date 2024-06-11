% Radar test at Cat Creek Summit 12/23/2022
% HP Marshall
% Radar by Adrian Tang et al, JPL and UCLA, as part of NASA IIP20
clear all; close all;

flow=15e9; % [Hz] start freq
fhigh=15.5e9; % [Hz] stop freq
BW=fhigh-flow; % [Hz] bandwidth
Tpl=67.55e-6; % [s] pulse length
Fs=122.88e6; % [Hz] sample rate
v=3.0e8; % [m/s] speed in air
N=2^15; % number of points in FFT
w=(0:N/2-1)/(N)*Fs; % frequencies sampled
d=0.5*w*Tpl/(BW)*v;

%% Sky calibration Variables
skyfile = 'P:\SnowDrones\Surveys\2024\2024-03-27_GrandMesa\Radar\Ku\Skycal_60m.csv';
crange=[3 20]; % Adrians initial values
thresh=2e9;
dsmooth=1;


%% Load a profile
D=load('P:\SnowDrones\Surveys\2024\2024-03-27_GrandMesa\Radar\Ku\data1.csv');
D=D'; % transpose due to collection in row vectors
D=D(:,1:length(D));
%D=D(:,1:390);
[nr,nc]=size(D);

%% Calibration scheme
S=load(skyfile); % load sky cal
[nrS,ncS] = size(S); % get size
% Crop the skycal if needed
if ncS > nr
    % Crop S to match the number of rows in D
    S = S(:, 1:nr);
    [nrS, ncS] = size(S);
elseif ncS < nr
    % Handle the case where S has fewer rows than D
    error('Array S has fewer rows (%d) than array D (%d).', nrS, nr);
end
average_cal = mean(S(10:nrS-10,:))'; % remove/start and end traces, to avoid weird stuff
PCAL = zeros(nr,nc); % we overwrite below, so changed this to just initialize (HPM)
TCAL = zeros(nr,nc); % changed to remove the original from the average (HPM)
%apply a sliding normalization point copy of the skycal
for calpoint = crange(1):crange(2) % range to calibrate over
    for i = 1:nc % loop traces
        fact = D(calpoint,i)./average_cal(calpoint); % normalize to calpoint
        PCAL(:,i) = D(:,i) - average_cal .* fact; % subtract average skycal trace, normalized
    end
    TCAL=TCAL+PCAL; % update matrix (sum for all calpoints)
end
TCAL=imgaussfilt(TCAL,dsmooth); % smooth in 2D
TCAL(TCAL<thresh)=0; % threshold amplitudes, set to zero if less than thresh


%% Plot the data
ix=1:nc; d2=d(1:nr);
figure(1);clf
% Plot the original data
subplot(2,1,1);
imagesc(ix,d2, D); colorbar;  % Optional, adds a color bar to indicate the color scale
title('Flight Transect', 'FontSize', 18);
subtitle('Original Radar Traces', 'FontSize', 14);
xlabel('Trace Number', 'FontSize', 16);
ylabel('Distance in Air [m]', 'FontSize', 16);
% caxis(cRange);
axis tight;
% Plot the skycal data
subplot(2,1,2);
imagesc(ix,d2,TCAL); colorbar
title('Flight Transect', 'FontSize', 18);
subtitle('Sky Calibration Applied', 'FontSize', 14);
xlabel('Trace Number', 'FontSize', 16); 
ylabel('Distance in Air [m]', 'FontSize', 16);





