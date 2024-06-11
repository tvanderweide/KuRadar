%% EMTS2 Pit Comparison
% Script for importing geotagged radar traces and selecting the traces
% closest to the pit
% v2 Adds a button to view the traces in a new window to save for the
% presentation
% Thomas Van Der Weide
% 5/20/2024

%% Initializing Variables
clear all; close all;
% Radar
flow=15e9; % [Hz] start freq
fhigh=15.5e9; % [Hz] stop freq
BW=fhigh-flow; % [Hz] bandwidth
Tpl=67.55e-6; % [s] pulse length
Fs=122.88e6; % [Hz] sample rate
v=3.0e8; % [m/s] speed in air
N=2^15; % number of points in FFT
w=(0:N/2-1)/(N)*Fs; % frequencies sampled
d=0.5*w*Tpl/(BW)*v;

% Global
global combinedData;
global pitLoc;
global pitLocFold;
global radarHeightInd;
global waveformDist;
global vSnow;
global overallDensityAverage;
global overallTemperatureAverage;
global pitHeight;
global pitName;
global pitDate;
global pitNameDate;
global radarFN;
global nr;
global annotationHandles;
global pitFile;
global pitFold;
global S;


%% Load the pit location csv
% Prompt user to specify the pit Location CSV
[pitLocFile, pitLocFold] = uigetfile({'*.csv';'*.CSV';'*.*'}, 'Select a Pit Data File','P:\SnowDrones\Surveys\2024\');
slashBeforePitDataIndex = strfind(pitLocFold, 'PIT_DATA\') - 1;
radarFoldStart = pitLocFold(1:slashBeforePitDataIndex);
% Load the selected File
if isequal(pitLocFile,0)
    disp('User selected Cancel');
else
    disp(['User will open ', fullfile(pitLocFold, pitLocFile)]);
    pitLocFilePath = fullfile(pitLocFold, pitLocFile);
end
Loc_opts = detectImportOptions(pitLocFilePath);
pitLoc = readtable(pitLocFilePath, Loc_opts);


%% Load the geotagged radar traces and skycal file
%  Prompt user to specify the pit location CSV
if isequal(radarFoldStart,0)
    [file, fold] = uigetfile({'*_geotagged.mat'}, 'Select a Radar Trace CSV File','P:\SnowDrones\Surveys\2024\');
else
    [file, fold] = uigetfile({'*_geotagged.mat'}, 'Select a Radar Trace CSV File',radarFoldStart);
end

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
traceMatrix = cell2mat(combinedData.RadarTrace.'); % transpose due to collection in row vectors
[nr,nc]=size(traceMatrix);
waveformDist = d(1:nr);


%% Calculate distances between Pit and radar traces
% Constants for conversion
earthRadius = 6371000; % Earth radius in meters
latConv = earthRadius * pi/180; % Conversion factor for latitude to meters
lonConv = @(lat) cosd(lat) * latConv; % Conversion factor for longitude to meters depending on latitude

% Loop through each pit location
for i = 1:height(pitLoc)
    % Calculate conversion factor for longitude based on pitLoc's latitude
    currentLonConv = lonConv(pitLoc.Lat(i));
    
    % Calculate Cartesian distances
    deltaX = (combinedData.longitudedeg - pitLoc.Lon(i)) * currentLonConv;
    deltaY = (combinedData.latitudedeg - pitLoc.Lat(i)) * latConv;

    % Euclidean distance
    distances = sqrt(deltaX.^2 + deltaY.^2);

    % Add distances as a new column to combinedData with dynamic name based on pitID
    columnName = sprintf('minDistanceTo%s', pitLoc.PitID{i});
    combinedData.(columnName) = distances;
end

%%
interactive_trace_plot()
function interactive_trace_plot()
    % Define variables
    global combinedData;
    global pitLoc;
    global subplotIndex;
    global distCalc;
    global radarHeightInd;
    global waveformDist;
    global vSnow;
    global radRes;
    global overallDensityAverage;
    global overallTemperatureAverage;
    global ixLowMin;
    global ixLowMax;
    global currentPitName;
    global pitHeight;
    global pitNameDate;
    global pitName;
    global pitDate;
    global radarFN;
    global annotationHandles;
    global radarDistance;
    global pitFile;
    global pitFold;
    global traceMatrixSubset;
    global sfile;
    global dsmooth;
    global thresh;
    global crange;
    global mrho;
    
    distCalc = 15; % Distance in meters from locations to select points
    subplotIndex = 1;
    
    D = cell2mat(combinedData.RadarTrace.');
    sfile='P:\SnowDrones\Surveys\2024\2024-03-27_GrandMesa\Radar\Ku\skycal_60m.csv'; % Skycal
    mrho=0; % for air
    dsmooth=1; thresh=50; crange=[3 20]; % inputs to procJPLradar.m
    [Z,Zp,x,waveformDist]=procJPLradar2(D,sfile,mrho,dsmooth,thresh,crange); % load and process data with skycal
    traceMatrix = Zp;
    pitNames = combinedData.Properties.VariableNames(contains(combinedData.Properties.VariableNames, 'minDistanceTo'));
    
    % Setup the figure and maximize it
    figure(2); clf;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    % Prepare six subplots
    % show Radar Traces in the first subplot
    subplot(2,3,1:3);
    imagesc(x,waveformDist,traceMatrix, [0 4e10]); colorbar;  % Optional, adds a color bar to indicate the color scale
    title(sprintf('Radar Trace Data (%s)',radarFN), 'FontSize', 18);
%     title('1-30-24 Flight1 Transect', 'FontSize', 18);
    subtitle(sprintf('Pits Located within %s [m]', num2str(distCalc)), 'FontSize', 14);
%     ylabel('Sample Index');
    ylabel('Distance in Air [m]', 'FontSize', 16);
    xlabel('Trace Number', 'FontSize', 16);
    axis tight;  % Fit the axes tightly around the data
    hold on;
    for i = 1:length(pitNames)
        simpleName = strrep(pitNames{i}, 'minDistanceTo', '');
        % Find indices where distance is less than 10 meters
        closeIndices = find(combinedData.(pitNames{i}) <= distCalc);
%         subsetCloseIndices = closeIndices(1:3:end);
         % Plot vertical lines at these indices
         if length(closeIndices) > 10
            xline(closeIndices(1:3), 'Color',[1 1 1 0.7], 'Label', simpleName, 'LabelOrientation', 'horizontal', 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom', 'FontSize', 14); % Plot a colored vertical line
            xline(closeIndices(end-2:end), 'Color',[1 0 0 0.7], 'Label', simpleName, 'LabelOrientation', 'horizontal', 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'bottom', 'FontSize', 14); % Plot a colored vertical line
         end
    hold off;
    end
    
    handles.ax4 = subplot(2, 3, 4);
    handles.ax5 = subplot(2, 3, 5);
    handles.ax6 = subplot(2, 3, 6);
    handles.newMax = 0.5 * max(traceMatrix(:));
    handles.legendData = [];
    set(gcf, 'UserData', handles);
    
    % Create UI buttons
    %Swap between selected Pits
    uicontrol('Style', 'pushbutton', 'String', '<', 'Units', 'normalized', ...
              'Position', [0.17 0.48 0.05 0.05], 'Callback', @(~,~) updateSubplots(-1, pitNames, combinedData, pitLoc));
    uicontrol('Style', 'pushbutton', 'String', '>', 'Units', 'normalized', ...
              'Position', [0.24 0.48 0.05 0.05], 'Callback', @(~,~) updateSubplots(1, pitNames, combinedData, pitLoc));
    annotation('textbox', [0.02 0.02 0.15 0.05], 'String', sprintf('Max Distance from Pit Location: %.2f', distCalc), ...
                   'EdgeColor', 'none', 'FontSize', 10, 'HorizontalAlignment', 'right');
    % Slider to adjust color scale
    uicontrol('Style', 'slider', 'Min',0,'Max',1,'Value',0.5, 'Units', 'normalized', ...
              'Position', [0.4 0.05 0.2 0.03], 'Callback', @adjustColorScale);
    txtMaxColorLimit = uicontrol('Style', 'text', 'Units', 'normalized', ...
              'Position', [0.45 0.025 0.2 0.02], 'String', sprintf('Max Color: %.2e', handles.newMax));
    % Subset radar Trace Data
    annotation('textbox', [0.45 0.5 0.05 0.03], 'String', "-", ...
                   'EdgeColor', 'none', 'FontSize', 10, 'HorizontalAlignment', 'right');
    hMin = uicontrol('Style', 'edit', 'Units', 'normalized', ...
              'Position', [0.43 0.49 0.05 0.03], 'String', '0');
    hMax = uicontrol('Style', 'edit', 'Units', 'normalized', ...
              'Position', [0.485 0.49 0.05 0.03], 'String', '0');
    uicontrol('Style', 'pushbutton', 'String', 'Update', 'Units', 'normalized', ...
              'Position', [0.54 0.49 0.05 0.03], 'Callback', @updateXLines);
    
    % Create text input and button for displayWindow
    uicontrol('Style', 'text', 'String', 'Dsmooth', 'Units', 'normalized', ...
              'Position', [0.01 0.44 0.03 0.02]);
    hDsmooth = uicontrol('Style', 'edit', 'Units', 'normalized', ...
                         'Position', [0.04 0.44 0.03 0.02], 'String', '1');
    uicontrol('Style', 'text', 'String', 'Thresh', 'Units', 'normalized', ...
              'Position', [0.01 0.48 0.03 0.02]);
    hThresh = uicontrol('Style', 'edit', 'Units', 'normalized', ...
                        'Position', [0.04 0.48 0.03 0.02], 'String', '50');
    uicontrol('Style', 'text', 'String', 'Crange', 'Units', 'normalized', ...
              'Position', [0.01 0.52 0.03 0.02]);
    hCrange1 = uicontrol('Style', 'edit', 'Units', 'normalized', ...
                         'Position', [0.04 0.52 0.04 0.02], 'String', '3');
    hCrange2 = uicontrol('Style', 'edit', 'Units', 'normalized', ...
                         'Position', [0.08 0.52 0.04 0.02], 'String', '20');
    uicontrol('Style', 'pushbutton', 'String', 'Update', 'Units', 'normalized', ...
              'Position', [0.01 0.40 0.07 0.03], ...
              'Callback',  @(src, event) updateRadarSubset(src, event, hDsmooth, hThresh, hCrange1, hCrange2));
    hDisplayLow = uicontrol('Style', 'edit', 'Units', 'normalized', ...
                               'Position', [0.622 0.09 0.05 0.03], 'String', '20');
    hDisplayHigh = uicontrol('Style', 'edit', 'Units', 'normalized', ...
                               'Position', [0.622 0.06 0.05 0.03], 'String', string(length(waveformDist)));
    uicontrol('Style', 'pushbutton', 'String', 'Update Display', 'Units', 'normalized', ...
              'Position', [0.622 0.03 0.05 0.03], 'Callback', @updateRadarDepthPlot);
    % Pit Data buttons
    uicontrol('Style', 'pushbutton', 'String', 'Change Pit', 'Units', 'normalized', ...
              'Position', [0.92 0.9 0.05 0.03], 'Callback', @readPitData);
    annotationHandles.date = annotation('textbox', [0.83 0.85 0.15 0.05], 'String', sprintf('Pit Date: %s', pitDate), ...
                   'EdgeColor', 'none', 'FontSize', 10, 'HorizontalAlignment', 'right');
    annotationHandles.id = annotation('textbox', [0.83 0.82 0.15 0.05], 'String', sprintf('Pit ID: %s', pitName), ...
                   'EdgeColor', 'none', 'FontSize', 10, 'HorizontalAlignment', 'right');
    annotationHandles.height = annotation('textbox', [0.83 0.79 0.15 0.05], 'String', sprintf('Pit Height [cm]: %d', round(pitHeight)), ...
                   'EdgeColor', 'none', 'FontSize', 10, 'HorizontalAlignment', 'right');
    annotationHandles.density = annotation('textbox', [0.83 0.76 0.15 0.05], 'String', sprintf('Pit Density [kg/m^3]: %d', round(overallDensityAverage)), ...
           'EdgeColor', 'none', 'FontSize', 10, 'HorizontalAlignment', 'right');
    annotationHandles.temp = annotation('textbox', [0.83 0.73 0.15 0.05], 'String', sprintf('Pit Temp [C]: %.1f', overallTemperatureAverage), ...
                   'EdgeColor', 'none', 'FontSize', 10, 'HorizontalAlignment', 'right');
    annotationHandles.velocity = annotation('textbox', [0.83 0.70 0.15 0.05], 'String', sprintf('V_S [m/s]: %.2e', vSnow), ...
                   'EdgeColor', 'none', 'FontSize', 10, 'HorizontalAlignment', 'right');
    annotationHandles.res = annotation('textbox', [0.83 0.67 0.15 0.05], 'String', sprintf('Resolution [cm]: %.1f', radRes), ...
                   'EdgeColor', 'none', 'FontSize', 10, 'HorizontalAlignment', 'right');
    annotationHandles.height2 = annotation('textbox', [0.83 0.18 0.15 0.05], 'String', sprintf('Pit Height [cm]: %d', round(pitHeight)), ...
                   'EdgeColor', 'none', 'FontSize', 10, 'HorizontalAlignment', 'right');
    annotationHandles.hYCoordText = annotation('textbox', [0.83 0.15 0.15 0.05], 'String', 'Radar Depth [cm]: 000', ...
                   'EdgeColor', 'none', 'FontSize', 10, 'HorizontalAlignment', 'right');
    uicontrol('Style', 'pushbutton', 'String', 'Save Files', 'Units', 'normalized', ...
              'Position', [0.92 0.13 0.05 0.03], 'Callback', @saveFiles);
    % Create text input and button for displayWindow
    hLineLow = uicontrol('Style', 'edit', 'Units', 'normalized', ...
                               'Position', [0.92 0.09 0.05 0.03], 'String', '0');
    hLineHigh = uicontrol('Style', 'edit', 'Units', 'normalized', ...
                               'Position', [0.92 0.06 0.05 0.03], 'String', '0');
    uicontrol('Style', 'pushbutton', 'String', 'Update Lines', 'Units', 'normalized', ...
              'Position', [0.92 0.03 0.05 0.03], 'Callback', @(src, event) updateLinesFromInput(hLineLow, hLineHigh, handles.ax6));
    updateSubplots(0, pitNames, combinedData, pitLoc);  % Initial plot
    % Create new figure to save
    uicontrol('Style', 'pushbutton', 'String', 'New Figure', 'Units', 'normalized', ...
          'Position', [0.01 0.01 0.05 0.03], 'Callback', @createNewFigure);

    function adjustColorScale(src, ~)
        newMax = src.Value * max(traceMatrixSubset(:));
        handles.newMax = newMax;
        caxis(handles.ax5, [0 newMax]); % Adjust color limits based on slider value and maximum intensity
        set(txtMaxColorLimit, 'String', sprintf('Max Color: %.2e', newMax)); % Update the text box to show new max color limit
        set(gcf, 'UserData', handles);
    end
    
    function updateSubplots(direction, pitNames, combinedData, pitLoc)
        % Update subplotIndex
        subplotIndex = subplotIndex + direction;
        if subplotIndex < 1
            subplotIndex = length(pitNames);
        elseif subplotIndex > length(pitNames)
            subplotIndex = 1;
        end
        
        % Clear previous plots in the second row of subplots
        for k = 4:6
            subplot(2, 3, k);
            cla;
        end
        
        % Get current pit name and find close indices
        currentPitName = pitNames{subplotIndex};
        closeIndices = find(combinedData.(currentPitName) < distCalc);
        assignin('base', 'closeIndices', closeIndices);
        % Fourth subplot: Plot Lat/Lon of the pit and all closeIndices
        subplot(2, 3, 4);
        pitID = strrep(currentPitName, 'minDistanceTo', '');
        pitRow = find(strcmp(pitLoc.PitID, pitID));
        nCloseIndices = length(closeIndices);
        if ~isempty(pitRow)
            pitLat = pitLoc.Lat(pitRow);
            pitLon = pitLoc.Lon(pitRow);
            plot(pitLon, pitLat, 'gs', 'MarkerSize', 10, 'LineWidth', 2);
            hold on;
            plot(combinedData.longitudedeg(closeIndices), combinedData.latitudedeg(closeIndices), 'bo');
            hold off;
            title(sprintf('Location: %s (n = %d)', pitID, nCloseIndices));
            xlabel('Longitude');
            ylabel('Latitude');
            legend({'Pit Location', 'Radar Trace Location'}, 'Location', 'best');
        else
            title('Pit Location Not Found');
        end
    
        % Plot the Radar Traces that are closest to the Pit
        dsmooth=1; thresh=50; crange=[3 20];
        [Z,Zp,x,waveformDist]=procJPLradar2(D,sfile,mrho,dsmooth,thresh,crange); % load and process data with skycal
        traceMatrixSubset = Zp(:,closeIndices);
        subplot(2, 3, 5);
        imagesc(traceMatrixSubset,[0 2e10]); colorbar;  % Optional, adds a color bar to indicate the color scale
        title('Subset of Radar Trace Data');
        xlabel('Trace Number'); ylabel('Sample Index');
        axis tight;  % Fit the axes tightly around the data
                
        % Sixth subplot can be updated with other desired data
        subplot(2, 3, 6);
        title('Radar Depth vs Pit Depth');
    end

    function updateRadarSubset(src, event, hDsmooth, hThresh, hCrange1, hCrange2)
        % Retrieve values from UI
        dsmooth = str2double(get(hDsmooth, 'String'));
        thresh = str2double(get(hThresh, 'String'));
        crange = [str2double(get(hCrange1, 'String')), str2double(get(hCrange2, 'String'))];
    
        % Plot the Radar Traces that are closest to the Pit
        [Z,Zp,x,waveformDist]=procJPLradar2(D,sfile,mrho,dsmooth,thresh,crange); % load and process data with skycal
        traceMatrixSubset = Zp(:,closeIndices);
        subplot(2, 3, 5);
        imagesc(traceMatrixSubset,[0 handles.newMax]); colorbar;  % Optional, adds a color bar to indicate the color scale
        title('Subset of Radar Trace Data');
        xlabel('Trace Number'); ylabel('Sample Index');
        axis tight;  % Fit the axes tightly around the data
%         % Overlay radarHeightInd
%         hold on;  % Keep the current plot
%         plot(1:length(closeIndices), radarHeightInd(closeIndices), 'w', 'LineWidth', 2);  % Plot radarHeightInd as a white line
%         hold off;  % Release the plot hold
    end

    function updateXLines(~, ~)
        % Fetch values from text boxes
        ixLowMin = str2double(get(hMin, 'String'));
        ixLowMax = str2double(get(hMax, 'String'));
        % Generate an array of positions to place x-lines
        ixLowArr = ixLowMin:2:ixLowMax;
        hold(handles.ax5, 'on');  % Ensure that the image is not erased
        delete(findall(handles.ax5, 'Type', 'ConstantLine'));  % Remove existing lines
        % Plot new x-lines for each position in ixLowArr
        arrayfun(@(x) xline(handles.ax5, x, 'r'), ixLowArr);
        hold(handles.ax5, 'off');

        % Update subplots
        updateRadarLocPlot()
        updateRadarDepthPlot()
    end

    function updateRadarLocPlot(~,~)
        % Load variables
        ixLowMin = str2double(get(hMin, 'String'));
        ixLowMax = str2double(get(hMax, 'String'));
        ixLowArr = ixLowMin:ixLowMax;
        
        % Update fourth subplot so selected indices display as blue
        subplot(2, 3, 4);
        pitID = strrep(currentPitName, 'minDistanceTo', '');
        pitRow = find(strcmp(pitLoc.PitID, pitID));
        nCloseIndices = length(closeIndices);
        if ~isempty(pitRow)
            pitLat = pitLoc.Lat(pitRow);
            pitLon = pitLoc.Lon(pitRow);
            plot(pitLon, pitLat, 'gs', 'MarkerSize', 10, 'LineWidth', 2);
            hold on;
            plot(combinedData.longitudedeg(closeIndices), combinedData.latitudedeg(closeIndices), 'bo');
            plot(combinedData.longitudedeg(closeIndices(ixLowArr)), combinedData.latitudedeg(closeIndices(ixLowArr)), 'ro');
            hold off;
            title(sprintf('Location: %s (n = %d)', pitID, nCloseIndices));
            xlabel('Longitude');
            ylabel('Latitude');
            legend({'Pit Location', 'Radar Trace Location', 'Selected traces'}, 'Location', 'best');
        else
            title('Pit Location Not Found');
        end
    end

    function updateRadarDepthPlot(~,~)
        % Load variables
        displayLow = str2double(get(hDisplayLow, 'String'));
        displayHigh = str2double(get(hDisplayHigh, 'String'));
        if isnan(displayLow) || displayLow <= 0
            displayLow = 10;  % Reset to default if input is invalid
            set(hDisplayLow, 'String', '10');  % Reset the text in the input box
        end
        if isnan(displayHigh) || displayHigh > length(waveformDist)
            displayHigh = length(waveformDist);  % Reset to default if input is invalid
            set(hDisplayHigh, 'String', string(length(waveformDist)));  % Reset the text in the input box
        end
        ixLowMin = str2double(get(hMin, 'String'));
        ixLowMax = str2double(get(hMax, 'String'));
        ixLowArr = ixLowMin:ixLowMax;

        % Function to update the sixth plot, assume all necessary variables are available
        cla(handles.ax6);
        firstLocation = traceMatrixSubset(:,ixLowArr);
        M1 = median(firstLocation, 2);  % Median of selected traces
        plot(handles.ax6, M1, waveformDist, 'r-', 'LineWidth', 3);  % Redraw the plot
        set(handles.ax6, 'YDir', 'reverse', 'LineWidth', 2, 'FontSize', 14);
        ylim(handles.ax6, [waveformDist(displayLow) waveformDist(displayHigh)]);
        title(handles.ax6, 'Median of Selected Radar Traces');
        xlabel(handles.ax6, 'Power spectral density');
        ylabel(handles.ax6, ['Distance in snow [m] using $\rho=' sprintf('%.0f', mrho) '[kg/m^3]$'], 'Interpreter', 'latex');
        % Set the ButtonDownFcn for the plot
        set(handles.ax6, 'ButtonDownFcn', @(src, event) selectYIndex(src, event, handles.ax6));
end


    function selectYIndex(src, event, ax)
        % Get the current point of the mouse click on the plot
        cp = get(ax, 'CurrentPoint');
        yCoord = cp(1, 2);
    
        % Store the y-indices as global variables
        persistent yIndMin yIndMax
        if isempty(yIndMin) || ~isempty(yIndMax)
            % Reset the indices if both are set
            yIndMin = [];
            yIndMax = [];
            % Clear previous horizontal lines
            delete(findall(ax, 'Type', 'line', 'Tag', 'clickLine'));
        end
    
        if isempty(yIndMin)
            yIndMin = findClosestIndex(waveformDist, yCoord);
            % Draw horizontal line at the selected y-coordinate
            line(ax, get(ax, 'XLim'), [waveformDist(yIndMin) waveformDist(yIndMin)], 'Color', 'b', 'LineWidth', 2, 'Tag', 'clickLine');
            set(hLineLow, 'String', num2str(yIndMin));
            set(annotationHandles.hYCoordText, 'String', 'Radar Depth [cm]: 000');
        elseif isempty(yIndMax)
            yIndMax = findClosestIndex(waveformDist, yCoord);
            % Draw horizontal line at the selected y-coordinate
            line(ax, get(ax, 'XLim'), [waveformDist(yIndMax) waveformDist(yIndMax)], 'Color', 'b', 'LineWidth', 2, 'Tag', 'clickLine');
            set(hLineHigh, 'String', num2str(yIndMax));
            % Update the text box with the Radar Depth difference
            radarDistance = abs(waveformDist(yIndMin) - waveformDist(yIndMax))*100;
            set(annotationHandles.hYCoordText, 'String', ['Radar Depth [cm]: ', num2str(radarDistance, '%.1f')]);
        end
    end

    function updateLinesFromInput(hLineLow, hLineHigh, ax)
        % Get the values from the text inputsRadar Depth
        yIndMin = str2double(get(hLineLow, 'String'));
        yIndMax = str2double(get(hLineHigh, 'String'));
    
        % Clear previous horizontal lines
        delete(findall(ax, 'Type', 'line', 'Tag', 'clickLine'));
    
        % Draw horizontal lines at the selected y-coordinates
        line(ax, get(ax, 'XLim'), [waveformDist(yIndMin) waveformDist(yIndMin)], 'Color', 'b', 'LineWidth', 2, 'Tag', 'clickLine');
        line(ax, get(ax, 'XLim'), [waveformDist(yIndMax) waveformDist(yIndMax)], 'Color', 'b', 'LineWidth', 2, 'Tag', 'clickLine');
    
        % Update the text box with the Radar Depth difference
        radarDistance = abs(waveformDist(yIndMin) - waveformDist(yIndMax))*100;
        set(annotationHandles.hYCoordText, 'String', ['Radar Depth [cm]: ', num2str(radarDistance, '%.1f')]);
    end

    function idx = findClosestIndex(arr, value)
        % Find the index of the closest value in the array
        [~, idx] = min(abs(arr - value));
    end

    function saveFiles(~, ~)    
        % Define the file paths
%         outFigureBase = fullfile(pitFold, 'RadarDepthvsPit');
        outFigureBase = ['P:\SnowDrones\Surveys\2024\Processing\', pitNameDate];
        outFigureFullPath = [outFigureBase, '.png'];
        counter = 1;
        while exist(outFigureFullPath, 'file')
            outFigureFullPath = sprintf('%s_%d.png', outFigureBase, counter);
            counter = counter + 1;
        end
        
%         outFileBase = fullfile(pitFold, [pitName, '_', radarFN]);
        outFileBase = ['P:\SnowDrones\Surveys\2024\Processing\', [pitNameDate, '_', radarFN]];
        outFileFullMat = [outFileBase, '.mat'];
        counter = 1;  % Initialize a counter
        while exist(outFileFullMat, 'file')
            outFileFullMat = sprintf('%s_%d.mat', outFileBase, counter);
            counter = counter + 1;
        end
        dataPath = 'P:\SnowDrones\Surveys\2024\2024_datav3.csv';
        
        % Save the current figure (Figure 2) as a PNG file
        bigFig = figure(2);
        saveas(bigFig, outFigureFullPath);
        disp(['Figure saved to: ', outFigureFullPath]);
    
        % Check if the CSV file exists
        if exist(dataPath, 'file')
            % File exists, append the data
            fid = fopen(dataPath, 'a');  % Open the file for appending
        else
            % File does not exist, create file and write the data
            fid = fopen(dataPath, 'w');  % Open the file for writing
            fprintf(fid, 'PitDate,PitID,Pit_Height,Radar_Distance, RadarFN,Density,yIndMin,yIndMax,fileName,selectedMin,selectedMax\n');  % Write header
        end
        % Write or append the data; handle NaN and numerical data properly
        yIndMin = str2double(get(hLineLow, 'String'));
        yIndMax = str2double(get(hLineHigh, 'String'));
        if isnan(pitHeight)
            fprintf(fid, '%s,%s,%d,%.1f,%s,%.0f,%d,%d,%s\n', pitDate, pitName, '0', radarDistance, radarFN,mrho,yIndMin,yIndMax,[pitNameDate, '_', radarFN],ixLowMin,ixLowMax);  % Leave pitHeight blank if NaN
        else
            fprintf(fid, '%s,%s,%d,%.1f,%s,%.0f,%d,%d,%s\n', pitDate, pitName, round(pitHeight), radarDistance, radarFN,mrho,yIndMin,yIndMax,[pitNameDate, '_', radarFN],ixLowMin,ixLowMax);  % Write data
        end
         fclose(fid);
        % Provide feedback to the user
        disp(['Data saved to ' dataPath]);
        % Save the selected traces
        save(outFileFullMat, 'traceMatrixSubset');
        disp(['Data saved to: ', outFileFullMat]);
        disp(['Finished Writing Data and Saving Figure']);
    end
    
    % Load Pit Data and calculate velocity in snow
    function readPitData(~,~)
        global pitLocFold nr;
        flow=15e9; % [Hz] start freq
        fhigh=15.5e9; % [Hz] stop freq
        BW=fhigh-flow; % [Hz] bandwidth
        Tpl=67.55e-6; % [s] pulse length
        Fs=122.88e6; % [Hz] sample rate
        v=3.0e8; % [m/s] speed in air
        N=2^15; % number of points in FFT
        w=(0:N/2-1)/(N)*Fs; % frequencies sampled
    
        % Prompt user to specify the pit xlsx 
        if isequal(pitLocFold,0)
            [pitFile, pitFold] = uigetfile({'*.xlsx';'*.XLSX';'*.*'}, 'Select a Pit Data File','P:\SnowDrones\Surveys\2024\');
        else
            [pitFile, pitFold] = uigetfile({'*.xlsx';'*.XLSX';'*.*'}, 'Select a Pit Data File',pitLocFold);
        end
        
        if isequal(pitFile,0)
            disp('User selected Cancel');
        else
            pitNameDate = extractBefore(pitFile,'.');
            underscoreIndex = strfind(pitNameDate, '_');
            if ~isempty(underscoreIndex)
                pitName = extractBefore(pitNameDate, underscoreIndex(1));
                pitDate = extractAfter(pitNameDate, underscoreIndex(1));
            end
        
            % Load selected file
            disp(['User will open ', fullfile(pitFold, pitFile)]);
            pitFilePath = fullfile(pitFold, pitFile);
            % Import the data
            % Get information about the sheets in the workbook
            [~, sheets] = xlsfinfo(pitFilePath);
            % Check if 'FRONT' sheet exists or if there's only one sheet
            if any(strcmp(sheets, 'FRONT'))
                sheetToUse = 'FRONT';
            elseif length(sheets) == 1
                sheetToUse = sheets{1};  % Use the only sheet available
            else
                error('No suitable sheet found, and multiple sheets are present.');
            end
            opts = detectImportOptions(pitFilePath, 'Sheet', sheetToUse);
            pit_data = readtable(pitFilePath, opts);
            pit_data(:, 11:end) = [];
            pit_data(:,2) = [];
            pit_data(:,2) = [];
            % Find the first row where "HeightAboveGround_1" equals 0
            zeroRow = find(pit_data.HeightAboveGround_1 == 0, 1, 'first');
            % If zero is not found, find the first NaN
            if isempty(zeroRow)
                zeroRow = find(isnan(pit_data.HeightAboveGround_1), 1, 'first');
            end
            % If a zero or NaN is found, delete all rows from that row onwards
            if ~isempty(zeroRow)
                pit_data(zeroRow:end, :) = [];
            end
            % Handle potential non-numeric entries safely
            numericA = str2double(pit_data.DensityProfileA);
            numericB = str2double(pit_data.DensityProfileB);
            numericExtra = str2double(pit_data.ExtraDensity);
            numericT = str2double(pit_data.T);
            % Check and ensure no size mismatch occurs
            if length(numericA) == height(pit_data)
                pit_data.DensityProfileA = numericA;
            end
            if length(numericB) == height(pit_data)
                pit_data.DensityProfileB = numericB;
            end
            if length(numericExtra) == height(pit_data)
                pit_data.ExtraDensity = numericExtra;
            end
            if length(numericT) == height(pit_data)
                pit_data.T = numericT;
            end
            averageDensityValues = nanmean([pit_data.DensityProfileA, pit_data.DensityProfileB, pit_data.ExtraDensity], 2);
            overallDensityAverage = nanmean(averageDensityValues);
            mrho=overallDensityAverage;
            overallTemperatureAverage = nanmean([pit_data.T]);
            pitHeight = pit_data.HeightAboveGround_1(1);
            e_s = e_snowdry(overallDensityAverage, flow, overallTemperatureAverage);
            e = real(e_s);
            vSnow =3.0e8/(sqrt(e)); % [m/s] speed in snow
            d2=0.5*w*Tpl/(BW)*vSnow;
%             waveformDist = d2(1:nr);
            [Z,Zp,x,waveformDist]=procJPLradar2(D,sfile,mrho); % load and process data with skycal
            radRes = abs(waveformDist(2) - waveformDist(1))*100;
            % Update the text in the figure
            set(annotationHandles.date, 'String', sprintf('Pit Date: %s', pitDate));
            set(annotationHandles.id, 'String', sprintf('Pit ID: %s', pitName));
            set(annotationHandles.height, 'String', sprintf('Pit Height [cm]: %d', round(pitHeight)));
            set(annotationHandles.height2, 'String', sprintf('Pit Height [cm]: %d', round(pitHeight)));
            set(annotationHandles.density, 'String', sprintf('Pit Density [kg/m^3]: %d', round(overallDensityAverage)));
            set(annotationHandles.temp, 'String', sprintf('Pit Temp [C]: %.1f', overallTemperatureAverage));
            set(annotationHandles.velocity, 'String', sprintf('V_S [m/s]: %.2e', vSnow));
            set(annotationHandles.res, 'String', sprintf('Resolution [cm]: %.1f', radRes));

            
        end
    end

    function createNewFigure(~, ~)
        figure(3); clf;
        set(gcf, 'Units', 'normalized', 'Position', [0, 0, 1, 0.5]);
        ax2_new = subplot(1,13,2:6);
        ax3_new = subplot(1,13,8:12);

% %         % Plot Lat/Lon of the pit and all closeIndices
% %         % Fetch values from text boxes
        ixLowMin = str2double(get(hMin, 'String'));
        ixLowMax = str2double(get(hMax, 'String'));
% %         ixLowArr = ixLowMin:ixLowMax;
% %         pitID = strrep(currentPitName, 'minDistanceTo', '');
% %         pitRow = find(strcmp(pitLoc.PitID, pitID));
% %         nCloseIndices = length(closeIndices);
% %         if ~isempty(pitRow)
% %             pitLat = pitLoc.Lat(pitRow);
% %             pitLon = pitLoc.Lon(pitRow);
% %             plot(ax1_new, pitLon, pitLat, 'gs', 'MarkerSize', 10, 'LineWidth', 2);
% %             hold(ax1_new, 'on');
% %             plot(ax1_new, combinedData.longitudedeg(closeIndices), combinedData.latitudedeg(closeIndices), 'bo');
% %             plot(ax1_new, combinedData.longitudedeg(closeIndices(ixLowArr)), combinedData.latitudedeg(closeIndices(ixLowArr)), 'ro');
% %             hold(ax1_new, 'off');
% %             title(ax1_new, sprintf('Location: %s (n = %d)', pitID, nCloseIndices));
% %             xlabel(ax1_new, 'Longitude');
% %             ylabel(ax1_new, 'Latitude');
% %             legend(ax1_new, {'Pit Location', 'Radar Trace Locations', 'Selected Traces'}, 'Location', 'best');
% %         end
        
        % Plot the Subset of radar traces
        dsmooth = str2double(get(hDsmooth, 'String'));
        thresh = str2double(get(hThresh, 'String'));
        crange = [str2double(get(hCrange1, 'String')), str2double(get(hCrange2, 'String'))];
        mrho = overallDensityAverage;
        % Plot the Radar Traces that are closest to the Pit
        [Z,Zp,x,waveformDist]=procJPLradar2(D,sfile,mrho,dsmooth,thresh,crange); % load and process data with skycal
        traceMatrixSubset = Zp(:,closeIndices);
        imagesc(ax2_new,traceMatrixSubset,[0 handles.newMax]); colorbar;  % Optional, adds a color bar to indicate the color scale
        title(ax2_new, 'Radar Trace Data');
        xlabel(ax2_new, 'Trace Number'); ylabel(ax2_new, 'Sample Index');
%         axis tight;

        % Load variables for depth plot
        yIndMin = str2double(get(hLineLow, 'String'));
        yIndMax = str2double(get(hLineHigh, 'String'));
        snowDistance = waveformDist(yIndMin);
        dist = waveformDist - snowDistance;
        pitHeight_m = pitHeight / 100;
        [~, nearestIndex] = min(abs(dist - pitHeight_m));

        displayLow = yIndMin-20;
        displayHigh = yIndMax+20;
        if displayHigh > length(waveformDist)
            displayHigh = length(waveformDist);
        end

        % Generate an array of positions to place x-lines
        ixLowArr = ixLowMin:2:ixLowMax;
        hold(ax2_new, 'on');  % Ensure that the image is not erased
        delete(findall(ax2_new, 'Type', 'ConstantLine'));  % Remove existing lines
        % Plot new x-lines for each position in ixLowArr
        arrayfun(@(x) xline(ax2_new, x, 'r'), ixLowArr);
        hold(ax2_new, 'off');
        ylim(ax2_new, [yIndMin-20 yIndMax+20]);
%         xlim(ax2_new, [x(ixLowMin-50) x(ixLowMax+50)])


        % Function to update the radar median plot
        ixLowArr = ixLowMin:ixLowMax;
        cla(ax3_new);
        firstLocation = traceMatrixSubset(:,ixLowArr);
        M1 = median(firstLocation, 2);  % Median of selected traces
        hMedian = plot(ax3_new, M1, dist, 'k-', 'LineWidth', 3);  % Redraw the plot
        set(ax3_new, 'YDir', 'reverse', 'LineWidth', 2, 'FontSize', 14);
        ylim(ax3_new, [dist(displayLow) dist(displayHigh)]);
        title(ax3_new, 'Selected Radar Traces and Measured Snow Depth');
        xlabel(ax3_new, 'Power spectral density');
        ylabel(ax3_new, ['Distance in snow [m] using $\rho=' sprintf('%.0f', mrho) '[kg/m^3]$'], 'Interpreter', 'latex');
        % Drawing a horizontal lines
        hline_surface = line(ax3_new, xlim(ax3_new), [dist(yIndMin) dist(yIndMin)], 'Color', 'blue', 'LineStyle', '-', 'LineWidth', 2);
        hline_ground = line(ax3_new, xlim(ax3_new), [dist(yIndMax) dist(yIndMax)], 'Color', 'red', 'LineStyle', '-', 'LineWidth', 2);
        hline_pitDepth = line(ax3_new, xlim(ax3_new), [dist(nearestIndex) dist(nearestIndex)], 'Color', 'green', 'LineStyle', '--', 'LineWidth', 2);
        legend(ax3_new, [hMedian, hline_surface,hline_ground, hline_pitDepth], {'Median of Radar Traces', 'Radar Surface Return', 'Radar Ground Return' , 'Measured Pit Depth'}, 'Location', 'northeast');

        
        % Replicate annotations on the left side of the new figure
        annotation('textbox', [0.001 0.8 0.15 0.05], 'String', get(annotationHandles.date, 'String'), ...
                   'EdgeColor', 'none', 'FontSize', 16, 'HorizontalAlignment', 'left');
        annotation('textbox', [0.001 0.74 0.15 0.05], 'String', get(annotationHandles.id, 'String'), ...
                   'EdgeColor', 'none', 'FontSize', 16, 'HorizontalAlignment', 'left');
        annotation('textbox', [0.001 0.68 0.15 0.05], 'String', get(annotationHandles.height, 'String'), ...
                   'EdgeColor', 'none', 'FontSize', 16, 'HorizontalAlignment', 'left');
        annotation('textbox', [0.001 0.62 0.15 0.05], 'String', get(annotationHandles.density, 'String'), ...
                   'EdgeColor', 'none', 'FontSize', 16, 'HorizontalAlignment', 'left');
        annotation('textbox', [0.001 0.56 0.15 0.05], 'String', get(annotationHandles.temp, 'String'), ...
                   'EdgeColor', 'none', 'FontSize', 16, 'HorizontalAlignment', 'left');
        annotation('textbox', [0.001 0.49 0.15 0.05], 'String', get(annotationHandles.velocity, 'String'), ...
                   'EdgeColor', 'none', 'FontSize', 16, 'HorizontalAlignment', 'left');
        annotation('textbox', [0.001 0.43 0.15 0.05], 'String', get(annotationHandles.res, 'String'), ...
                   'EdgeColor', 'none', 'FontSize', 16, 'HorizontalAlignment', 'left');

        annotation('textbox', [0.83 0.18 0.15 0.05], 'String', get(annotationHandles.height2, 'String'), ...
                   'EdgeColor', 'none', 'FontSize', 16, 'HorizontalAlignment', 'right');
        radarDistance = abs(waveformDist(yIndMin) - waveformDist(yIndMax))*100;
        annotation('textbox', [0.83 0.12 0.15 0.05], 'String', ['Radar Depth [cm]: ', num2str(radarDistance, '%.1f')], ...
                   'EdgeColor', 'none', 'FontSize', 16, 'HorizontalAlignment', 'right');
    end

end


