%% EMTS2 Geotagging
% Script for importing radar GPS position file and geotagging radar traces
% Thomas Van Der Weide
% 4/3/2024

%% Variables to change
clear all; close all;
global TCAL;
global D;
global radarPosInterp;
global outFileFull;
global outFileFullMat;
global outFigureFullPath;
global intervalms;
intervalms = 0;

% Prompt user to specify a file location for saving
[gpsFile, gpsFold] = uigetfile({'*.POS';'*.txt';'*.*'}, 'Select a GPS Pos File','P:\SnowDrones\Surveys\2024\');
if isequal(gpsFile,0)
    disp('User selected Cancel');
%     gpsFold = 'P:\SnowDrones\Surveys\2024\2024-03-25_GrandMesa\GPSDATA\ReachRadar_raw_20240325174921_UBX\';
%     gpsFile = 'ReachRadar_raw_20240325174921.pos';
else
    disp(['User will open ', fullfile(gpsFold, gpsFile)]);
    % Save data to the specified file, for example:
    % save(fullfile(path, file), 'dataVariable');
end
S=load('P:\SnowDrones\Surveys\2024\2024-03-27_GrandMesa\Radar\Ku\skycal_60m.csv'); % Skycal

[file, fold] = uigetfile({'*.csv'}, 'Select a Radar Trace CSV File','P:\SnowDrones\Surveys\2024\');
if isequal(file,0)
    disp('User selected Cancel');
%     fold = 'P:\SnowDrones\Surveys\2024\2024-03-25_GrandMesa\Radar\Ku\'; % Radar trace folder
%     file = 'data1.csv'; % Radar Traces
else
    disp(['User will open ', fullfile(fold, file)]);
    % Save data to the specified file, for example:
    % save(fullfile(path, file), 'dataVariable');
end
outFileFull =  [fold, strrep(file, '.csv', '_geotagged.csv')];
outFileFullMat =  [fold, strrep(file, '.csv', '_geotagged.mat')];
outFigureFullPath = [fold,strrep(file, '.csv', '_Figure.png')];

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


%% load a profile
D=load(append(fold, file));
D=D'; % transpose due to collection in row vectors
%D=D(:,1:390); % Crop the recorded distance
[nr,nc]=size(D);


%% Calibration scheme
[nrS, ncS] = size(S);
% Crop the skycal if needed
if ncS > nr
    % Crop S to match the number of rows in D
    S = S(:, 1:nr);
    [nrS, ncS] = size(S);
elseif ncS < nr
    % Handle the case where S has fewer rows than D
    error('Array S has fewer rows (%d) than array D (%d).', nrS, nr);
end
%average all the traces in the sky cal excluding the start and end ones
%because weird things always happen at the start and end of a dataset
average_cal = mean(S(10:ncS-10,:))';
S=S'; % transpose due to collection in row vectors

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

%% Import the radar GPS data and resample
gpsPOS = append(gpsFold,gpsFile);

opts = delimitedTextImportOptions("NumVariables", 15);
% Specify range and delimiter
opts.DataLines = [11, Inf];
opts.Delimiter = " ";
% Specify column names and types
opts.VariableNames = ["Date", "UTC", "latitudedeg", "longitudedeg", "heightm", "Q", "ns", "sdnm", "sdem", "sdum", "sdnem", "sdeum", "sdunm", "ages", "ratio"];
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
dateAndTime = datetime(strcat(char(radarPos.Date), " " ,char(radarPos.UTC)), 'InputFormat', 'yyyy/MM/dd HH:mm:ss.SSS');
radarPos.("DateTime") = dateAndTime;
radarPos.DateTime.Format = 'yyyy/MM/dd HH:mm:ss.SSS';

% Resample GPS data into 54ms intervals to match the radar sample rate
if ncS == 801
    intervalms = 54;
    new_time = radarPos.DateTime(1):milliseconds(intervalms):radarPos.DateTime(end);
    new_time = new_time';
    radarPosInterp = table(new_time);
elseif ncS == 401
    intervalms = 52;
    new_time = radarPos.DateTime(1):milliseconds(intervalms):radarPos.DateTime(end);
    new_time = new_time';
    radarPosInterp = table(new_time);
elseif ncS == 251
    intervalms = 50;
    new_time = radarPos.DateTime(1):milliseconds(intervalms):radarPos.DateTime(end);
    new_time = new_time';
    radarPosInterp = table(new_time);
else
    print("Radar trace length is not known.");
    return;
end

% interp each column
for col = 3:5
    values_interpolated = interp1(radarPos.DateTime, radarPos{:, col}, new_time, 'linear');
    radarPosInterp = [radarPosInterp, array2table(values_interpolated, 'VariableNames', {radarPos.Properties.VariableNames{col}})];
end

%% Use just the radar traces
% % Pull the radar traces
% traceMatrix = cell2mat(combinedData.RadarTrace.');
% % Plot
% figure;
% imagesc(traceMatrix);
% colorbar;  % Optional, adds a color bar to indicate the color scale
% title('Visualization of Radar Trace Data');
% xlabel('Sample Index');
% ylabel('Trace Number');
% axis tight;  % Fit the axes tightly around the data
% colormap jet; 


%%
interactive_trace_plot()

function interactive_trace_plot()
    global TCAL;  % Ensure TCAL is recognized in this function
    global D;
    global radarPosInterp;
    global minIdx;
    global maxIdx;
    global min_Radartrace;
    global max_Radartrace;
    global RadarTraceCount;
    global GPSCount;
    global DifferenceCount;
    global updateFlag;    % enable/disable updating comparison plot
    global posCropped;
    global tcal_data;
    global trace_data;
    global outFileFull;
    global outFileFullMat;
    global outFigureFullPath;
    global intervalms;
    
    % Initialize index counters
    minIdx = 1;
    maxIdx = numel(radarPosInterp.new_time);
    min_Radartrace = 1;
    max_Radartrace = size(TCAL, 2);
    RadarTraceCount = max_Radartrace;
    GPSCount = size(radarPosInterp,1);
    updateFlag = false; % Initially, update flag is false
    DifferenceCount = GPSCount - max_Radartrace;
    
    % Setup the figure and maximize it
    figure(2); clf;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

    % Prepare four subplots; initially show TCAL in the first two subplots
    subplot(3,2,1);
    imagesc(TCAL); colorbar;
    title('Full TCAL Data Visualization');
    xlabel('Trace Number');
    ylabel('Sample Number');

    subplot(3,2,2);
    imagesc(TCAL); colorbar;
    title('TCAL Data Visualization');
    xlabel('Trace Number');
    ylabel('Sample Number');

    % Third subplot remains for showing initial full radarPosInterp data
    subplot(3,2,3);
    plot(radarPosInterp.new_time, radarPosInterp.heightm, 'b-');
    title('Initial Full Height vs. Time');
    xlabel('Time');
    ylabel('Height (m)');
    grid on;
    datetick('x', 'keeplimits');

    % Fourth subplot for dynamic updates via sliders
    subplot(3,2,4);
    h4 = plot(radarPosInterp.new_time, radarPosInterp.heightm, 'b-'); % Initial plot
    title('Dynamic Height vs. Time');
    xlabel('Time');
    ylabel('Height (m)');
    grid on;

    % Bottom row subplot for dynamic updates
    subplot(3,2,[5,6]);
    h56 = plot(radarPosInterp.new_time, radarPosInterp.heightm, 'b-'); % Initial plot
    title('TCAL vs. GPS Time');
    xlabel('Time');
    ylabel('Height (m)');
    legend('gpsAlt','radarHeight');
    grid on;

    % Create UI controls for Radar Trace interaction
    RadartextBox = uicontrol('Style', 'edit', 'String', 'Radar Trace Selection', 'Position', [20 825 150 20]);
    minRadarBox = uicontrol('Style', 'edit', 'String', '1', 'Position', [20 800 150 20]);
    uicontrol('Style', 'pushbutton', 'String', '-', 'Position', [20, 775, 70, 20], 'Callback', @decrement_Radarmin);
    uicontrol('Style', 'pushbutton', 'String', '+', 'Position', [100, 775, 70, 20], 'Callback', @increment_Radarmin);
    maxRadarBox = uicontrol('Style', 'edit', 'String', string(length(TCAL)), 'Position', [20 750 150 20]);
    uicontrol('Style', 'pushbutton', 'String', '-', 'Position', [20, 725, 70, 20], 'Callback', @decrement_Radarmax);
    uicontrol('Style', 'pushbutton', 'String', '+', 'Position', [100, 725, 70, 20], 'Callback', @increment_Radarmax);
    uicontrol('Style', 'slider', 'Min', min_Radartrace, 'Max', max_Radartrace, 'Value', min_Radartrace, ...
              'Position', [20, 700, 150, 20], 'Callback', @update_Radar_plot_slider, 'Tag', 'minRadarSlider');
    uicontrol('Style', 'slider', 'Min', min_Radartrace, 'Max', max_Radartrace, 'Value', max_Radartrace, ...
              'Position', [20, 675, 150, 20],'Callback', @update_Radar_plot_slider, 'Tag', 'maxRadarSlider');
    uicontrol('Style', 'pushbutton', 'String', 'Update Radar Plot', 'Position', [20 650 150 20], 'Callback', @update_Radar_plot);
    RadarTraceCountDisplay = uicontrol('Style', 'edit', 'String', max_Radartrace, 'Position', [1650 800 150 20]);
    uicontrol('Style', 'edit', 'String', 'Radar Trace Count:', 'Position', [1650 825 150 20]);
    uicontrol('Style', 'edit', 'String', ['Intervalms: ', num2str(intervalms)], 'Position', [1650 900 150 20]);

    % Define the callback function to update the cropped TCAL subplot
    function update_Radar_plot(~, ~)
        min_Radartrace = str2double(get(minRadarBox, 'String'));  % Get the string from minBox and convert to double
        max_Radartrace = str2double(get(maxRadarBox, 'String'));  % Get the string from maxBox and convert to double
        RadarTraceCount = max_Radartrace - min_Radartrace + 1;
        DifferenceCount = GPSCount - RadarTraceCount;
        set(RadarTraceCountDisplay, 'String', num2str(RadarTraceCount));
        set(DifferenceCountDisplay, 'String', num2str(DifferenceCount));
        % Compare values and change background color
        if str2double(get(GPSCountDisplay, 'String')) == str2double(get(RadarTraceCountDisplay, 'String'))
            set(GPSCountDisplay, 'BackgroundColor', 'green');
            set(RadarTraceCountDisplay, 'BackgroundColor', 'green');
            set(DifferenceCountDisplay, 'BackgroundColor', 'green');
        else
            set(GPSCountDisplay, 'BackgroundColor', 'white');
            set(RadarTraceCountDisplay, 'BackgroundColor', 'white');
            set(DifferenceCountDisplay, 'BackgroundColor', 'white');
        end
        
        % Check if the input values are within the valid range
        if isnan(min_Radartrace) || isnan(max_Radartrace) || min_Radartrace < 1 || max_Radartrace > size(TCAL, 2) || min_Radartrace > max_Radartrace
            errordlg('Please enter valid trace numbers within the correct range.', 'Invalid Input', 'modal');
            return;
        end

        % Update only the second subplot
        subplot(3,2,2); 
        hold off;
        imagesc(TCAL(:, max(1,min_Radartrace):min(size(TCAL,2),max_Radartrace))); colorbar;
        title(['TCAL Data Visualization: Columns ' num2str(min_Radartrace) ' to ' num2str(max_Radartrace)]);
        xlabel('Trace Number');
        ylabel('Sample Number');
        % Update h56 subplot
        if updateFlag
            plotComp();
        end
    end

    function update_Radar_plot_slider(~, ~)       
        % Get the slider values
        min_Radartrace = round(get(findobj('Tag', 'minRadarSlider'), 'Value'));
        max_Radartrace = round(get(findobj('Tag', 'maxRadarSlider'), 'Value'));
        
        % Update the edit boxes
        set(minRadarBox, 'String', num2str(min_Radartrace));
        set(maxRadarBox, 'String', num2str(max_Radartrace));
        
        % Update the radar plot
        update_Radar_plot();
    end
    
    % Callback functions for increment/decrement buttons for min_Radartrace
    function decrement_Radarmin(~, ~)
        if min_Radartrace > 1
            min_Radartrace = min_Radartrace - 1;
            minRadarBox.String = num2str(min_Radartrace);
            update_Radar_plot()
        end
    end

    function increment_Radarmin(~, ~)
        if min_Radartrace < max_Radartrace
            min_Radartrace = min_Radartrace + 1;
            minRadarBox.String = num2str(min_Radartrace);
            update_Radar_plot()
        end
    end

    % Callback functions for increment/decrement buttons for max_Radartrace
    function decrement_Radarmax(~, ~)
        if max_Radartrace > min_Radartrace
            max_Radartrace = max_Radartrace - 1;
            maxRadarBox.String = num2str(max_Radartrace);
            update_Radar_plot()
        end
    end

    function increment_Radarmax(~, ~)
        if max_Radartrace < size(TCAL, 2)
            max_Radartrace = max_Radartrace + 1;
            maxRadarBox.String = num2str(max_Radartrace);
            update_Radar_plot()
        end
    end


    % Create UI controls for GPS Height vs. Time interaction
    sliderWidth = 150; % Width of each slider
    sliderHeight = 20; % Height of each slider
    spaceBetween = 20; % Space between sliders
    uicontrol('Style', 'edit', 'String', 'GPS Time Range Adjustment', 'Position', [20, 525, 150, 20]);
    minGPSBox = uicontrol('Style', 'edit', 'String', char(radarPosInterp.new_time(minIdx)), 'Position', [20, 500, 150, 20]);
    uicontrol('Style', 'pushbutton', 'String', '-', 'Position', [20, 475, 70, 20], 'Callback', @decrement_min);
    uicontrol('Style', 'pushbutton', 'String', '+', 'Position', [100, 475, 70, 20], 'Callback', @increment_min);
    maxGPSBox = uicontrol('Style', 'edit', 'String', char(radarPosInterp.new_time(maxIdx)), 'Position', [20, 450, 150, 20]);
    uicontrol('Style', 'pushbutton', 'String', '-', 'Position', [20, 425, 70, 20], 'Callback', @decrement_max);
    uicontrol('Style', 'pushbutton', 'String', '+', 'Position', [100, 425, 70, 20], 'Callback', @increment_max);
    uicontrol('Style', 'slider', 'Min', minIdx, 'Max', maxIdx, 'Value', minIdx, ...
              'Position', [20, 400, sliderWidth, sliderHeight], ...
              'Callback', @update_GPS_plot, 'Tag', 'minSlider');
    uicontrol('Style', 'slider', 'Min', minIdx, 'Max', maxIdx, 'Value', maxIdx, ...
              'Position', [20, 375, sliderWidth, sliderHeight], ...
              'Callback', @update_GPS_plot, 'Tag', 'maxSlider');
    uicontrol('Style', 'edit', 'String', 'GPS Time Entries:', 'Position', [1650 525 150 20]);
    GPSCountDisplay = uicontrol('Style', 'edit', 'String', GPSCount, 'Position', [1650 500 150 20]);
    % First row: Checkbox and Label for minIdx
    minIdxCheckbox = uicontrol('Style', 'checkbox', 'Position', [20, 350, 20, 20], 'Value', 0);
    uicontrol('Style', 'text', 'String', 'minIdx', 'Position', [45, 350, 50, 20]);
    % Second row: Checkbox and Label for maxIdx
    maxIdxCheckbox = uicontrol('Style', 'checkbox', 'Position', [20, 320, 20, 20], 'Value', 0);
    uicontrol('Style', 'text', 'String', 'maxIdx', 'Position', [45, 320, 50, 20]);
    % Third row: Entry box for entering the value to adjust the index
    indexValueEntry = uicontrol('Style', 'edit', 'String', '0', 'Position', [20, 290, 100, 20]);
    % Fourth row: Button to update the GPS index
    uicontrol('Style', 'pushbutton', 'String', 'Update GPS Index', 'Position', [20, 260, 150, 20], ...
    'Callback', @updateGPSIndex);

    % Function to update GPS index based on selected checkbox and entry value
    function updateGPSIndex(~, ~)
        value = str2double(get(indexValueEntry, 'String')); % Get the numeric value from the entry box
        if isnan(value) % Check if the input is not a number
            errordlg('Please enter a valid numeric value.');
            return;
        end
    
        if get(minIdxCheckbox, 'Value') == 1 % Check if minIdx is selected
            minIdx = max(1, minIdx + value); % Update minIdx, ensure it doesn't go below 1
            set(minGPSBox, 'String', char(radarPosInterp.new_time(minIdx))); % Update display
        end
    
        if get(maxIdxCheckbox, 'Value') == 1 % Check if maxIdx is selected
            maxIdx = max(minIdx, maxIdx + value); % Update maxIdx, ensure it doesn't go below minIdx
            set(maxGPSBox, 'String', char(radarPosInterp.new_time(maxIdx))); % Update display
        end
    
        % Redraw the plot with new indices
        update_plot();
    end

    % Callback function to update the fourth subplot based on slider values
    function update_GPS_plot(~, ~)
        minSlider = findobj(gcf, 'Tag', 'minSlider');
        maxSlider = findobj(gcf, 'Tag', 'maxSlider');
        minIdx = round(get(minSlider, 'Value'));
        maxIdx = round(get(maxSlider, 'Value'));
        valid_idx = minIdx:maxIdx;
        GPSCount = length(radarPosInterp.new_time(valid_idx));
        DifferenceCount = GPSCount - RadarTraceCount;
        % Update text boxes
        set(GPSCountDisplay, 'String', num2str(GPSCount));
        set(minGPSBox, 'String', char(radarPosInterp.new_time(minIdx)));
        set(maxGPSBox, 'String', char(radarPosInterp.new_time(maxIdx)));
        set(DifferenceCountDisplay, 'String', num2str(DifferenceCount));
        % Compare values and change background color
        if str2double(get(GPSCountDisplay, 'String')) == str2double(get(RadarTraceCountDisplay, 'String'))
            set(GPSCountDisplay, 'BackgroundColor', 'green');
            set(RadarTraceCountDisplay, 'BackgroundColor', 'green');
            set(DifferenceCountDisplay, 'BackgroundColor', 'green');
        else
            set(GPSCountDisplay, 'BackgroundColor', 'white');
            set(RadarTraceCountDisplay, 'BackgroundColor', 'white');
            set(DifferenceCountDisplay, 'BackgroundColor', 'white');
        end
        
        % Filter data
        subplot(3,2,4);
        set(h4, 'XData', radarPosInterp.new_time(valid_idx), 'YData', radarPosInterp.heightm(valid_idx));
        title('GPS Height vs. Time for Selected Range');
        xlabel('Time');
        ylabel('Height (m)');
        grid on;
        datetick('x', 'keeplimits'); % Adjust the x-axis ticks

        if updateFlag
            plotComp();
        end
    end
    
    % Callback functions for increment/decrement buttons
    function decrement_min(~, ~)
        if minIdx > 1
            minIdx = minIdx - 1;
            set(minGPSBox, 'String', char(radarPosInterp.new_time(minIdx)));
        end
        update_plot()
    end

    function increment_min(~, ~)
        if minIdx < numel(radarPosInterp.new_time)
            minIdx = minIdx + 1;
            set(minGPSBox, 'String', char(radarPosInterp.new_time(minIdx)));
        end
        update_plot()
    end

    function decrement_max(~, ~)
        if maxIdx > 1
            maxIdx = maxIdx - 1;
            set(maxGPSBox, 'String', char(radarPosInterp.new_time(maxIdx)));
        end
        update_plot()
    end
    
    function increment_max(~, ~)
        if maxIdx < numel(radarPosInterp.new_time)
            maxIdx = maxIdx + 1;
            set(maxGPSBox, 'String', char(radarPosInterp.new_time(maxIdx)));
        end
        update_plot()
    end
    % Callback function to update the plot using global counters
    function update_plot(~, ~)
        % Update text boxes
        set(minGPSBox, 'String', char(radarPosInterp.new_time(minIdx)));
        set(maxGPSBox, 'String', char(radarPosInterp.new_time(maxIdx)));
        valid_idx = minIdx:maxIdx;
        GPSCount = length(radarPosInterp.new_time(valid_idx));
        DifferenceCount = GPSCount - RadarTraceCount;
        % Compare values and change background color
        set(GPSCountDisplay, 'String', num2str(GPSCount));
        set(DifferenceCountDisplay, 'String', num2str(DifferenceCount));
        if str2double(get(GPSCountDisplay, 'String')) == str2double(get(RadarTraceCountDisplay, 'String'))
            set(GPSCountDisplay, 'BackgroundColor', 'green');
            set(RadarTraceCountDisplay, 'BackgroundColor', 'green');
            set(DifferenceCountDisplay, 'BackgroundColor', 'green');
        else
            set(GPSCountDisplay, 'BackgroundColor', 'white');
            set(RadarTraceCountDisplay, 'BackgroundColor', 'white');
            set(DifferenceCountDisplay, 'BackgroundColor', 'white');
        end
        

        % Update plot
        subplot(3,2,4);
        set(h4, 'XData', radarPosInterp.new_time(valid_idx), 'YData', radarPosInterp.heightm(valid_idx));
        title('GPS Height vs. Time for Selected Range');
        xlabel('Time');
        ylabel('Height (m)');
        grid on;
        datetick('x', 'keeplimits'); % Adjust the x-axis ticks

        % Adjust y-axis limits to enhance visibility of minimum values
        min_height = min(radarPosInterp.heightm(valid_idx));
        max_height = max(radarPosInterp.heightm(valid_idx));
        ylim([min_height - 0.1*(max_height - min_height), max_height]); % Add some padding to the y-axis

        if updateFlag
            plotComp();
        end
    end

    % Create UI controls for GPS Height vs. Radar Height
    uicontrol('Style', 'edit', 'String', 'Difference between counts:', 'Position', [1650 225 150 20]);
    DifferenceCountDisplay = uicontrol('Style', 'edit', 'String', DifferenceCount, 'Position', [1650 200 150 20]);
    % Add checkbox UI control
    updateCheckbox = uicontrol('Style', 'checkbox', 'String', 'Update Comparison Plot', 'Value', 0, 'Position', [20, 200, 20, 20], 'Callback', @toggleUpdateFlag);
    uicontrol('Style', 'text', 'String', 'Update Comparison Plot', 'Position', [40, 197, 100, 50], 'HorizontalAlignment', 'left');
%     uicontrol('Style', 'pushbutton', 'String', 'Combine GPS and Radar Traces', 'Position', [1650 175 150 20], 'Callback', @geotag);
    uicontrol('Style', 'pushbutton', 'String', 'Combine GPS and Radar Traces', 'Position', [1650 175 150 20], 'Callback', @(src, evnt) geotagCallback());

    % Callback function for the checkbox to toggle update flag
    function toggleUpdateFlag(source, ~)
        updateFlag = logical(source.Value);
    end

    function geotagCallback(~, ~)
        try
            combinedData = geotag();  % Call geotag with necessary arguments
            assignin('base', 'combinedData', combinedData);
        catch ME
            errordlg(ME.message, 'Error Combining Data');  % Display an error dialog with the message
        end
%         % Prompt the user to select a save location
%         [outFile, outPath] = uiputfile('*.csv', 'Save As');
%         % Check if the user canceled the operation
%         if isequal(outFile,0) || isequal(outPath,0)
%             disp('User canceled the operation');
%             return;
%         end
%         % Construct the full file path
%         outFileFull = fullfile(outPath, outFile);
        disp(['Writing Data', "..."]);
        % Write the table to a CSV file
        writetable(combinedData, outFileFull);
        disp(['Data saved to: ', outFileFull]);
        % Write the table to a MAT file
        save(outFileFullMat, 'combinedData');
        disp(['Data saved to: ', outFileFullMat]);
        % Save the Figure
        bigFig = figure(2);
        saveas(bigFig, outFigureFullPath);
        disp(['Figure saved to: ', outFigureFullPath]);
        disp(['Finished Writing Data', ""]);
    end

    function plotComp(~, ~)
        % Pull the appropriate GPS Data
        valid_idx = minIdx:maxIdx;
        posCropped=radarPosInterp(valid_idx,:);
        x = 1:size(posCropped.heightm,1);
        gpsHeight = posCropped.heightm;

        % Extract the radar height above ground
        valid_idxR = min_Radartrace:max_Radartrace; % Index range for TCAL
        tcal_data = TCAL(:, valid_idxR);
        trace_data = D(:, valid_idxR);
        [M,I] = max(tcal_data);
        radarHeightInd = movmean(I,40);
        
        % Scale the radar index to match GPS height
        gpsHeightDiff = max(gpsHeight) - min(gpsHeight);
        maxV = max(radarHeightInd);
        minV = min(radarHeightInd);
        radarHeight = ((radarHeightInd - minV) / (maxV - minV)) * gpsHeightDiff;
        radarAlt = radarHeight + min(posCropped.heightm) - 2;

        % Ensure sizes match by padding with zeros if needed
        if size(tcal_data, 2) < size(gpsHeight, 1)
            tcal_data = [tcal_data, zeros(size(tcal_data, 1), size(gpsHeight, 1) - size(tcal_data, 2))];
        elseif size(tcal_data, 2) > size(gpsHeight, 1)
            gpsHeight = [gpsHeight, zeros(size(gpsHeight, 1), size(tcal_data, 2) - size(gpsHeight, 1))];
        end
        
        % Plot
        subplot(3,2,[5,6]);
        plot(radarAlt, 'r');
        hold on;
        plot(x,posCropped.heightm, 'b');
        hold off;
        xlabel('Time');
        ylabel('Height above ground (m)');
        legend('radarHeight','GPS Height');
        title('Alignment of Radar Height and GPS Height');

    end

    function combinedData = geotag(~, ~)
        assignin('base', 'posCropped', posCropped);
        assignin('base', 'trace_data', trace_data);
        
        % Ensure that the number of columns in tcal_data matches the number of rows in posCropped
        if size(posCropped, 1) ~= size(trace_data, 2)
            error('Mismatch in the number of rows in posCropped and columns in tcal_data.');
        end
    
        % Convert each column of trace_data into a cell array where each cell contains one column vector
        trace_cells = mat2cell(trace_data, size(trace_data, 1), ones(1, size(trace_data, 2)));
        assignin('base', 'trace_cells', trace_cells);
        % Initialize combinedData from posCropped
        combinedData = posCropped;
        % Add this cell array as a new column in posCropped
        % Check if a specific column name is already used and adjust if necessary
        columnName = 'RadarTrace';
        if ismember(columnName, combinedData.Properties.VariableNames)
            error('The column name "%s" already exists in the table. Choose a different name.', columnName);
        end
        combinedData.(columnName) = trace_cells';
        
        % Check correct data was written to the file
        traceMatrix = cell2mat(combinedData.RadarTrace.');
        % Plot
        figure(3);clf;
        imagesc(traceMatrix);
        colorbar;  % Optional, adds a color bar to indicate the color scale
        title('Visualization of Radar Trace Data');
        xlabel('Sample Index');
        ylabel('Trace Number');
        axis tight;  % Fit the axes tightly around the data
        colormap jet; 

    end
end




