% procJPLradar
% HPM, improved with skycal from AT
% 4/16/2023
% process
% INPUT: dfile = csv file from radar
%        skyfile = csv file from radar pointed at sky
%        dsmooth = distance to smooth over
%        thresh = threshold, pixels with amplitude<thresh are set to 0
%        crange = calibration point range to use for skycal normalization 
% OUTPUT: Z = original radar matrix
%        Zp = processed radar matrix (skycal, smoothing, threshold)
%         x = trace numbers
%         dist = distance in air

function [Z,Zp,x,dist]=procJPLradar(dfile,skyfile,mrho,dsmooth,thresh,crange)
if nargin<6
    crange=[3 20]; % adrians initial values
end
if nargin<5
    thresh=1e9;
end
if nargin<4
    dsmooth=1;
end
if nargin<3
    mrho=0; % use distance in air
end

% define distance scale from radar parameters
flow=15e9; % [Hz] start freq
fhigh=15.996e9; % [Hz] stop freq
BW=fhigh-flow; % [Hz] bandwidth
Tpl=67.55e-6; % [s] pulse length
Fs=122.88e6; % [Hz] sample rate
v=3.0e8; % [m/s] speed in air
N=2^14; % number of points in FFT
w=(0:N/2-1)/(N)*Fs; % frequencies sampled
d=0.5*w*Tpl/(BW)*v; % distance in air
D=load(dfile); % load radar data file
Z=D'; % transpose due to collection in row vectors
[nr,nc]=size(Z); % get size
x=1:nc; dist=d(1:nr); % make vectors for plot
es=e_snowdry(mrho,15.5e9,-5); % get diel const for given density (mrho=0 for air)
vsnow=v./sqrt(real(es)); % snow velocity
dist=d(1:nr)*vsnow/v; % distance in snow

if nargin>1 % if skycal is given, use Adrian's approach to normalize and subtract
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
            fact = Z(calpoint,i)./average_cal(calpoint); % normalize to calpoint
            PCAL(:,i) = Z(:,i) - average_cal .* fact; % subtract average skycal trace, normalized
        end
        TCAL=TCAL+PCAL; % update matrix (sum for all calpoints)
    end
    TCAL=imgaussfilt(TCAL,dsmooth); % smooth in 2D
    TCAL(TCAL<thresh)=0; % threshold amplitudes, set to zero if less than thresh
end
Zp=TCAL; % output processed radar result
% uncomment below to plot
%h=imagesc(x,dist,Z,zrange); colorbar
