clear;
% Plot a comparison over the same Pit at two flight heights
% Thomas Van Der Weide and HP Marshall
% Radar by Adrian Tang et al, JPL and UCLA, as part of NASA IIP20

flow=15e9; % [Hz] start freq
fhigh=15.5e9; % [Hz] stop freq
BW=fhigh-flow; % [Hz] bandwidth
Tpl=67.55e-6; % [s] pulse length
Fs=122.88e6; % [Hz] sample rate
v=3.0e8; % [m/s] speed in air
N=2^15; % number of points in FFT
w=(0:N/2-1)/(N)*Fs; % frequencies sampled
d=0.5*w*Tpl/(BW)*v;

%% Take the sky calibration
%D=load('C:\Users\thomasvanderweide\Documents\RedPitaya\EMTS-2\2023-04-05_MCSLOWER\skycal.csv');
S=load('P:\SnowDrones\Surveys\2024\2024-03-27_GrandMesa\Radar\Ku\skycal_30m.csv');
[nr nc] = size(S);

%average all the traces in the sky cal excluding the start and end ones
%because weird things always happen at the start and end of a dataset
average_cal = mean(S(10:nr-10,:))';
S=S'; % transpose due to collection in row vectors


%% load a profile
fold = 'P:\SnowDrones\Surveys\2024\2024-03-26_GrandMesa\Radar\Ku\';
file = 'data1.csv';
D=load(append(fold, file));
D=D'; % transpose due to collection in row vectors
%D=D(:,1:390); % Crop the recorded distance
[nr,nc]=size(D);

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

%% Plot the entire file
transectDist = 1700; % How far the UAV transect was in meters
ixTrace=1:nc; iy=1:nr; d2=d(1:nr);
ixDist = linspace(1,transectDist,nc);
figure(2);clf;
imagesc(ixDist,d2,TCAL,[0 2e10]); colorbar
title(append('Full Transect from: ', file));

%% Refine the Plot
ixTrace=1:nc; iy=1:nr; d2=d(1:nr);
ixDist = linspace(1,transectDist,nc);
g = figure(3);clf;
subplot(1,12,1:7)
imagesc(ixDist,d2,TCAL,[0 2e10]); colorbar
%imagesc(ix,d2,TCAL,[0 2e6]); colorbar
%imagesc(ix,d2,TCAL,[0 2e2]); colorbar
grid on;
hold on; 
% Include vertical lines indicating where high and low traces were pulled
pitName = "Pit WENO13";
ixLow = 1045;
ixHigh = 1022;
window = 2; 
ixLowArr = linspace(ixLow-window,ixLow+window,window*2+1);
ixHighArr = linspace(ixHigh-window,ixHigh+window,window*2+1);
vLowLine = xline(ixLowArr, 'Color',[1 0 0 0.9]);
vHighLine = xline(ixHighArr, 'Color',[0 1 0 0.9]);
axis([(ixLow - 100) (ixHigh + 100) 10 45]); %Automatically adjust plot based on area of interest

title('3-27-2024 Grand Mesa ' + pitName);
xlabel('Position [m]'); ylabel('distance from radar in air [m]');
set(gca,'FontSize',30,'LineWidth',3,'FontWeight','bold');
% text(14,1.45,'snow surface','fontsize',24,'fontweight','bold','color','w')
% text(5,1.45,'snow surface','fontsize',24,'fontweight','bold','color','w')
% text(4,2.15,'ground','fontsize',24,'fontweight','bold','color','w')
% text(15,2.37,'ground','fontsize',24,'fontweight','bold','color','w')
% text(400,2.4,'low ground return - off nadir?','fontsize',14,'fontweight','bold','color','w')
% text(1600,1,'up-down','fontsize',18,'fontweight','bold','color','w')
% text(500,1,'~20 meter profile','fontsize',18,'fontweight','bold','color','w')

% Median of Traces with labels
% Find the correct index
lowerMin = ixLow-window;
lowerMax = ixLow+window;
[ valDiff, ixLowMin ] = min( abs( ixDist-lowerMin ) );
[ valDiff, ixLowMax ] = min( abs( ixDist-lowerMax ) );
upperMin = ixHigh-window;
upperMax = ixHigh+window;
[ valDiff, ixHighMin ] = min( abs( ixDist-upperMin ) );
[ valDiff, ixHighMax ] = min( abs( ixDist-upperMax ) );

firstLocation = D(:,ixLowMin:ixLowMax);
secondLocation = D(:,ixHighMin:ixHighMax);
M1=median(firstLocation,2); % get median at first location
M2=median(secondLocation,2); % get median at second location
subplot(1,12,9:12)
plot(M1,d2,'r-','linewidth',3); hold on
plot(M2,d2,'g-','linewidth',3);
legend('Lower Alt','Higher Alt');
set(gca, 'ydir','reverse');
%axis([0 5e8 0 60]);
ylim([10 45]);
set(gca,'linewidth',2,'fontsize',14);
title('Median of traces');
xlabel('power spectral density');
%ylabel('distance from radar in air [m]')
set(gca,'FontSize',30,'LineWidth',3,'FontWeight','bold');
%text(0.5e9,1.2,'snow surface','fontsize',14,'fontweight','bold','color','g');
%text(3e9,1.65,'snow surface','fontsize',14,'fontweight','bold','color','g');
%text(0.5e9,2.3,'ground?','fontsize',14,'fontweight','bold','color','g');
%text(0.6e9,1.75,'snow surface','fontsize',14,'fontweight','bold','color','r');
%text(3e9,2.1,'ground','fontsize',14,'fontweight','bold','color','r');

% Zoomed in figure
% Find the correct index
[ valDiff, ixLowMaxDisp ] = max(abs(M1-0));
[ valDiff, ixHighMaxDisp] = max(abs(M2-0));
displayWindow = 20;
% Plot
h = figure(4);clf
yyaxis left
plot(M1,d2,'r-','linewidth',3);
ylim([d2(ixLowMaxDisp-displayWindow) d2(ixLowMaxDisp+displayWindow+15)]);
set(gca, 'ydir','reverse');
ylabel('distance from radar in air [m]');
yyaxis right
plot(M2,d2,'g-','linewidth',3);
ylim([d2(ixHighMaxDisp-displayWindow) d2(ixHighMaxDisp+displayWindow+15)]);
legend('Lower Alt Flight','Higher Alt Hover');
title('Comparison of traces over ' + pitName);
ax = gca;
ax.YAxis(1).Color = 'r';
ax.YAxis(2).Color = 'g';
set(gca, 'ydir','reverse');
ylabel('distance from radar in air [m]');


%% Save Files
% "Saving Files to: "
% name = split(file,'.');
% B = erase(pitName," ");
% save(fold + B + "_" + string(name{1}) + "_" + string(round(d2(ixLowMaxDisp))) +"m","firstLocation")
% save(fold + B + "_" + string(name{1}) + "_" + string(round(d2(ixHighMaxDisp))) +"m","secondLocation")
% exportgraphics(g,fold + B + "_" + string(name{1}) + "_data.png")
% exportgraphics(h,fold + B + "_" + string(name{1}) + "_comp.png")

%% Save the Array for processing in Python
%save("C:\Users\Endure1\Documents\SnowDrones\Surveys\2024-01-29_GrandMesa\Radar\TCALdata2.mat","TCAL")
