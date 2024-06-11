%% EMTS2 Display a scatter plot the pit data results from 2024
% Thomas Van Der Weide
% 5/28/2024

% Read the CSV file
data = readtable('P:\SnowDrones\Surveys\2024\2024_Pit_Data.csv');

% Extract the fields
pit_height = data.Pit_Height;
radar_distance = data.Radar_Distance;
location = data.Location;
n = length(pit_height);

% Get unique locations and assign colors
unique_locations = unique(location);
colors = lines(length(unique_locations)); % Generate distinct colors

% Create a scatter plot
figure; hold on;
for i = 1:length(unique_locations)
    loc = unique_locations{i};
    idx = strcmp(location, loc);
    scatter(radar_distance(idx), pit_height(idx), 36, 'MarkerFaceColor', colors(i,:), 'DisplayName', loc);
end

% Fit a linear model for all data points
p = polyfit(radar_distance, pit_height, 1);
yfit = polyval(p, radar_distance);

% Plot the best fit line for all data points
plot(radar_distance, yfit, 'k-', 'LineWidth', 2, 'DisplayName', 'Best Fit Line');

% Add one-to-one line
minVal = min([radar_distance; pit_height]);
maxVal = max([radar_distance; pit_height]);
plot([minVal maxVal], [minVal maxVal], '--', 'Color', [0.5 0.5 0.5], 'DisplayName', '1:1 Line');


% Calculate R^2 value for all data points
yresid = pit_height - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(pit_height) - 1) * var(pit_height);
rsq = 1 - SSresid / SStotal;
rmse = sqrt(mean(yresid.^2));
% Display stats
xlims = xlim;
ylims = ylim;
text(xlims(1) + 0.05 * range(xlims), ylims(2) - 0.11 * range(ylims), ...
    {sprintf('N = %d', n),...
    sprintf('R^2 = %.2f', rsq),...
    sprintf('RMSE [cm] = %.2f', rmse)}, 'Color', 'k', 'FontSize', 12);

% Add labels and legend
xlabel('Radar Distance');
ylabel('Pit Height');
title('Pit Height vs. Radar Distance');
legend('Location', 'Best', 'Interpreter', 'none'); % Ensure underscores are not interpreted as subscripts
grid on; hold off;

