%% EMTS2 Geotagging
% Script for importing radar GPS position file and geotagging radar traces
% Thomas Van Der Weide
% 4/3/2024

%% Variables to change
clear all;
dateStr = "03-25-2024 ";
gpsFold = 'P:\SnowDrones\Surveys\2024\2024-03-25_GrandMesa\GPSDATA\ReachRadar_raw_20240325174921_UBX\';
gpsFile = 'ReachRadar_raw_20240325174921.pos';
S=load('P:\SnowDrones\Surveys\2024\2024-03-27_GrandMesa\Radar\Ku\skycal_60m.csv'); % Skycal
fold = 'P:\SnowDrones\Surveys\2024\2024-03-25_GrandMesa\Radar\Ku\'; % Radar trace folder
file = 'data1.csv'; % Radar Traces


%% Radar Variables
flow=15e9; % [Hz] start freq
fhigh=15.5e9; % [Hz] stop freq
BW=fhigh-flow; % [Hz] bandwidth
Tpl=67.55e-6; % [s] pulse length
Fs=122.88e6; % [Hz] sample rate
v=3.0e8; % [m/s] speed in air
N=2^15; % number of points in FFT
w=(0:N/2-1)/(N)*Fs; % frequencies sampled
d=0.5*w*Tpl/(BW)*v;


%% Import the radar GPS data and resample
gpsPOS = append(gpsFold,gpsFile);

opts = delimitedTextImportOptions("NumVariables", 15);
% Specify range and delimiter
opts.DataLines = [11, Inf];
opts.Delimiter = " ";
% Specify column names and types
opts.VariableNames = ["Date", "GPST", "latitudedeg", "longitudedeg", "heightm", "Q", "ns", "sdnm", "sdem", "sdum", "sdnem", "sdeum", "sdunm", "ages", "ratio"];
opts.VariableTypes = ["categorical", "string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";

% Specify variable properties
opts = setvaropts(opts, "Date", "EmptyFieldRule", "auto");

% Import the data
radarPos = readtable(gpsPOS, opts);
% Format the Date
dateAndTime = datetime(dateStr + radarPos.GPST, 'InputFormat', 'MM-dd-yyyy HH:mm:ss.SSS');
radarPos.("DateTime") = dateAndTime;
radarPos.DateTime.Format = 'MM-dd-yyyy HH:mm:ss.SSS';

% Resample GPS data into 54ms intervals to match the radar sample rate
new_time = radarPos.DateTime(1):milliseconds(54):radarPos.DateTime(end);
new_time = new_time';
radarPosInterp = table(new_time);

% interp each column
for col = 3:5
    values_interpolated = interp1(radarPos.DateTime, radarPos{:, col}, new_time, 'linear');
    radarPosInterp = [radarPosInterp, array2table(values_interpolated, 'VariableNames', {radarPos.Properties.VariableNames{col}})];
end


%% Take the sky calibration
[nr nc] = size(S);
%average all the traces in the sky cal excluding the start and end ones
%because weird things always happen at the start and end of a dataset
average_cal = mean(S(10:nr-10,:))';
S=S'; % transpose due to collection in row vectors

%% load a profile
D=load(append(fold, file));
D=D'; % transpose due to collection in row vectors
%D=D(:,1:390); % Crop the recorded distance
[nr,nc]=size(D);

%% Plot the entire radar trace file
transectDist = 1700; % How far the UAV transect was in meters
ixTrace=1:nc; iy=1:nr; d2=d(1:nr);
ixDist = linspace(1,transectDist,nc);
figure(1);clf;
imagesc(ixDist,d2,D,[0 2e9]); colorbar
title(append('Full Transect from: ', file));

%% Calibration scheme
%marker to prealign cal 
%D(37,:)=10e9;   %used to determine where the cal should be considered
PCAL = D;
TCAL = D;

%apply a sliding normalization point copy of the skycal
for calpoint = 3:20
  for i = 1:nc
  fact = D(calpoint,i)./average_cal(calpoint);
  PCAL(:,i) = D(:,i) - average_cal .* fact;
  end
TCAL=TCAL+PCAL;
end

%smoothing if you want it
TCAL=imgaussfilt(TCAL,0.7); 

%threshold the sliding cal
threshold = 5e9; %tunable threshold
map = TCAL<threshold;
TCAL(map)=0;

%% Sync the radar and GPS x-axis
figure(2);clf;
posCropped=radarPosInterp(33260:end-20805,:);
x = 1:size(posCropped.heightm,1);

% Extract the radar height above ground
[M,I] = max(TCAL);
radarHeightInd = movmean(I,40);
% Scale the radar index to match GPS height
gpsHeightDiff = max(posCropped.heightm) - min(posCropped.heightm);
maxV = max(radarHeightInd);
minV = min(radarHeightInd);
radarHeight = ((radarHeightInd - minV) / (maxV - minV)) * gpsHeightDiff;
radarAlt = radarHeight + min(posCropped.heightm) - 2;
% Plot
figure(3);clf;
plot(radarAlt);
hold on;
plot(x,posCropped.heightm);
legend('radarHeight','gpsAlt');


%% Geotag the radar traces
arr = table2array(posCropped);
test = D';
geotagged = [arr, test];

