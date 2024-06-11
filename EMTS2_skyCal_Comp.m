%% EMTS2 SkyCal Plot
% Script used to create a figure in the presentation showing original vs skyCal Radar Traces
% and a scatter plot of the pit data results from 2024
% Thomas Van Der Weide
% 6/3/2024

%% Initializing Variables
clear all; close all;

% Load radar files
sfile='P:\SnowDrones\Surveys\2024\2024-03-27_GrandMesa\Radar\Ku\skycal_60m.csv'; % Skycal

%  Prompt user to specify the pit location CSV
[file, fold] = uigetfile({'*_geotagged.mat'}, 'Select a Radar Trace CSV File','P:\SnowDrones\Surveys\2024\');

if isequal(file,0)
    disp('User selected Cancel');
%     fold = 'P:\SnowDrones\Surveys\2024\2024-03-25_GrandMesa\Radar\Ku\'; % Radar trace folder
%     file = 'data1.csv'; % Radar Traces
else
    disp(['User will open ', fullfile(fold, file)]);
    radarFN = extractBefore(file, "_");
    % Save data to the specified file, for example:
    % save(fullfile(path, file), 'dataVariable');
end

% load radar profile
load(append(fold, file));
D = cell2mat(combinedData.RadarTrace.');

%%
% Plotting params
dsmooth = 1;
thresh = 2e8;
crange = [3, 20];
mrho = 0;

figure(1);clf;
subplot(2,1,1);
[Z,Zp,x,waveformDist]=procJPLradar2(D,sfile,mrho,dsmooth,thresh,crange);
imagesc(x, waveformDist, Z, [0 2e9]); colorbar;  % Optional, adds a color bar to indicate the color scale
title('3-26-24 Flight1 Transect', 'FontSize', 18);
subtitle('Original Radar Traces', 'FontSize', 14);
xlabel('Trace Number', 'FontSize', 16);
ylabel('Distance in Air [m]', 'FontSize', 16);
% caxis(cRange);
axis tight;

subplot(2,1,2);
imagesc(x, waveformDist, Zp, [0 4e10]); colorbar;
title('3-26-24 Flight1 Transect', 'FontSize', 18);
subtitle('Sky Calibration Applied', 'FontSize', 14);
xlabel('Trace Number', 'FontSize', 16); 
ylabel('Distance in Air [m]', 'FontSize', 16);
% caxis(cRange);
axis tight;

% lowIdx = 2000;
% highIdx = 4000;
% closeIndicies = lowIdx:highIdx;
% traceMatrixSubsetZ = Z(:, closeIndicies);
% traceMatrixSubsetZp = Zp(:,closeIndicies);
% 
% displayLow = 235;
% displayHigh = 305;
% 
% subplot(2,2,2);
% imagesc(x(closeIndicies), waveformDist, traceMatrixSubsetZ); colorbar;
% title('3-26-24 Flight1 Transect Subset');
% subtitle('Sky Calibration Applied');
% xlabel('Trace Number'); 
% ylabel('Distance in Air [m]');
% ylim([waveformDist(displayLow) waveformDist(displayHigh)]);
% caxis(cRange);
% 
% subplot(2,2,4);
% imagesc(x(closeIndicies), waveformDist, traceMatrixSubsetZp); colorbar;
% ylim([waveformDist(displayLow) waveformDist(displayHigh)]);
% title('3-26-24 Flight1 Transect Subset');
% subtitle('Sky Calibration Applied');
% xlabel('Trace Number'); 
% ylabel('Distance in Air [m]');
% caxis(cRange);



