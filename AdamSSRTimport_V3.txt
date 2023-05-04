close all
clear
% clc

% Script to read multiple .lvm files into structure arrays

% Get a list of all .lvm files in the current directory
file_list = dir('*.lvm');

% Create a structure array to store all the file data
all_files_data = struct();

% Loop through each file and import the data into a structure array
for i = 1:length(file_list)
    % Get the file name without the extension
    [~, file_name] = fileparts(file_list(i).name);

    % Replace periods with underscores in the file name
    file_name = strrep(file_name, '.', '_');

    % Add a prefix to the file name to avoid issues with field names starting with numbers
    file_name = ['data_' file_name];

    % Initialize a temporary structure for the current file
    temp_struct = struct();

    % Open the file for reading
    fid = fopen(file_list(i).name, 'r');

    % Read header information
    for j = 1:21
        line = fgetl(fid);
        temp_struct.header{j} = line;
    end

    % Read column titles (Line 22)
    col_titles = strsplit(fgetl(fid), '\t');
    temp_struct.column_titles = col_titles;

    % Read the data (Line 23 onwards)
    data = textscan(fid, '%f %f %f %f', 'HeaderLines', 1, 'Delimiter', '\t', 'CollectOutput', 1);
    data = data{1};

    % Save data to the temporary structure
    temp_struct.X_Value = data(:, 1);
    temp_struct.Time = data(:, 2);
    temp_struct.Displacement = data(:, 3);
    temp_struct.Load = data(:, 4);

    % Close the file
    fclose(fid);

    % Convert Displacement from inches to mm and add as a new column
    temp_struct.Displacement = [temp_struct.Displacement, temp_struct.Displacement * 25.4];

    % Convert Load from lbf to N and add as a new column
    temp_struct.Load = [temp_struct.Load, temp_struct.Load * 4.44822];

    % Add new fields to the structure
    temp_struct.Gauge_Length = 27;
    temp_struct.Gauge_Width = 6.20;
    temp_struct.Engr_Stress = [];
    temp_struct.Engr_Strain = [];

%     % Prompt the user to input the Gauge_Thickness for each file
%     prompt = sprintf('Please input the sample thickness for the file %s:', file_name);
%     temp_struct.Gauge_Thickness = input(prompt);

    % Prompt the user to input the Gauge_Thickness for each file
    prompt = sprintf('Please input the sample thickness for the file %s:', file_name);
    temp_struct.Gauge_Thickness = 0;

    % Save the temporary structure to the all_files_data structure array
    all_files_data.(file_name) = temp_struct;
end

num_samples = length(file_list);

%Temporary Insert so I don't have to enter these fuckin' thickness. See
%above Prompt line to uncomment when done
%measurements during debugging%
all_files_data.data_5N_pure_10_h.Gauge_Thickness = 4.48;
all_files_data.data_5N_pure_10.Gauge_Thickness = 4.48;
all_files_data.data_5N_pure_30_h.Gauge_Thickness = 3.58;
all_files_data.data_5N_pure_30.Gauge_Thickness = 3.58;
all_files_data.data_5N_pure_60_h.Gauge_Thickness = 2.09;
all_files_data.data_5N_pure_60.Gauge_Thickness = 2.09;
%End Insert

% Add an additional sample '5N_pure_0' as a reference for undeformed sample state
reference_sample_name = '5N_pure_0';
all_files_data.data_5N_pure_0 = struct();

% Set Gauge_Length, Gauge_Width, and Gauge_Thickness for the reference sample
all_files_data.data_5N_pure_0.Gauge_Length = 27;
all_files_data.data_5N_pure_0.Gauge_Width = 6.20;
all_files_data.data_5N_pure_0.Gauge_Thickness = 5.01;

% Create empty fields for the reference sample
all_files_data.data_5N_pure_0.header = [];
all_files_data.data_5N_pure_0.column_titles = [];
all_files_data.data_5N_pure_0.X_Value = [];
all_files_data.data_5N_pure_0.Time = [];
all_files_data.data_5N_pure_0.Displacement = [];
all_files_data.data_5N_pure_0.Load = [];
all_files_data.data_5N_pure_0.Engr_Stress = [];
all_files_data.data_5N_pure_0.Engr_Strain = [];

% [...] Previous script parts go here

% Loop through each sample to offset Load and Displacement fields and calculate engineering stress and strain
sample_names = fieldnames(all_files_data);
for i = 1:length(sample_names)
    sample_name = sample_names{i};
    sample = all_files_data.(sample_name);

    % Check if the Load and Displacement fields are not empty
    if ~isempty(sample.Load) && ~isempty(sample.Displacement)
        % Calculate offset values for Load and Displacement fields
        initial_load = sample.Load(1, 2);
        initial_displacement = sample.Displacement(1, 2);

        sample.Load(:, 3) = sample.Load(:, 2) - initial_load;
        sample.Displacement(:, 3) = sample.Displacement(:, 2) - initial_displacement;

        % Replace negative values with the average of the cells in the row above and below
        for j = 2:length(sample.Load) - 1
            if sample.Load(j, 3) < 0
                sample.Load(j, 3) = (sample.Load(j - 1, 3) + sample.Load(j + 1, 3)) / 2;
            end
            if sample.Displacement(j, 3) < 0
                sample.Displacement(j, 3) = (sample.Displacement(j - 1, 3) + sample.Displacement(j + 1, 3)) / 2;
            end
        end

        % Calculate the cross-sectional area and store it in Gauge_xSection_Area_Initial
        sample.Gauge_xSection_Area_Initial = sample.Gauge_Thickness * sample.Gauge_Width;

        % Populate the Engr_Stress field by multiplying the offset Load (N) by the constant
        sample.Engr_Stress = sample.Load(:, 3) / sample.Gauge_xSection_Area_Initial;

        % Calculate and populate the Engr_Strain field
        sample.Engr_Strain = sample.Displacement(:, 3) / sample.Gauge_Length;
    else
        sample.Engr_Stress = [];
        sample.Engr_Strain = [];
    end

    % Update the sample in the all_files_data structure array
    all_files_data.(sample_name) = sample;
end

% Create a new figure for the interactive plot
figure('Name', 'Engineering Stress vs. Engineering Strain');
hold on;

% Set the X and Y labels for the plot
xlabel('Engr. Strain (mm/mm)');
ylabel('Engr. Stress (MPa)');

% Get the number of samples and their names
sample_names = fieldnames(all_files_data);

% Create arrays to store the plot handles and checkbox handles
plot_handles = gobjects(1, num_samples);

% Loop through each sample to create the plots and checkboxes
for i = 1:num_samples
    sample_name = sample_names{i};
    sample = all_files_data.(sample_name);

    % Create a plot for the current sample
    plot_handles(i) = plot(NaN, NaN, 'DisplayName', sample_name);
    set(plot_handles(i), 'Visible', 'off');

    % Update the plot data, accounting for different lengths of data
    set(plot_handles(i), 'XData', sample.Engr_Strain, 'YData', sample.Engr_Stress);
end

% Create a new figure for the checkboxes
checkbox_figure = figure('Name', 'Sample Selection', 'NumberTitle', 'off', ...
    'Position', [50, 50, 200, num_samples * 25 + 50]);

% Loop through each sample to create the checkboxes
for i = 1:num_samples
    sample_name = sample_names{i};
    display_name = strrep(sample_name, 'data_', '');

    % Create a checkbox for the current sample with the modified display name, horizontal alignment set to 'right'
    uicontrol('Style', 'checkbox', 'String', display_name, ...
        'HorizontalAlignment', 'left', 'FontSize', 8, 'Position', [10, checkbox_figure.Position(4) - 20 - i * 25, 180, 20], ...
        'Callback', {@toggle_plot, plot_handles(i)});
end



%**************************************************************************

% Load data
data = all_files_data.data_5N_pure_10;
stress = all_files_data.data_5N_pure_10.Engr_Stress;
strain = all_files_data.data_5N_pure_10.Engr_Strain;

% Perform k-means clustering
k = 3;
X = [strain, stress];
[idx, centroids] = kmeans(X, k);

% Identify cluster indices for each region
[~, compliance_cluster_idx] = min(centroids(:, 2));
[~, elastic_cluster_idx] = max(centroids(:, 1));
necking_cluster_idx = setdiff(1:k, [compliance_cluster_idx, elastic_cluster_idx]);

% Extract regions from the clustered data
compliance_region = X(idx == compliance_cluster_idx, :);
elastic_region = X(idx == elastic_cluster_idx, :);
necking_region = X(idx == necking_cluster_idx, :);

% Find the UTS and its index
[UTS, idx_UTS] = max(stress);

% Modify regions to only include points up to the UTS
compliance_region = compliance_region(compliance_region(:, 1) <= strain(idx_UTS), :);
elastic_region = elastic_region(elastic_region(:, 1) <= strain(idx_UTS), :);
necking_region = necking_region(necking_region(:, 1) <= strain(idx_UTS), :);

% Define plastic deformation and necking region (after UTS)
post_UTS_region = X(strain > strain(idx_UTS), :);

% % Plot the stress-strain curve with different colors for each region up to UTS
% figure;
% hold on;
% plot(compliance_region(:, 1), compliance_region(:, 2), 'b');
% plot(elastic_region(:, 1), elastic_region(:, 2), 'g');
% plot(necking_region(:, 1), necking_region(:, 2), 'r');
% plot(post_UTS_region(:, 1), post_UTS_region(:, 2), 'r'); % Plot post-UTS points as necking region (red)
% xlabel('Engineering Strain');
% ylabel('Engineering Stress');
% legend('Compliance', 'Elastic Deformation', 'Necking');

% Calculate the slope of the elastic region
elastic_slope = (elastic_region(end, 2) - elastic_region(1, 2)) / (elastic_region(end, 1) - elastic_region(1, 1));

% Print Young's modulus
fprintf("Young's Modulus: %.2f\n", elastic_slope);

% Plot the stress-strain curve with different colors for each region up to UTS
figure;
hold on;
plot(compliance_region(:, 1), compliance_region(:, 2), 'b');
plot(elastic_region(:, 1), elastic_region(:, 2), 'g');
plot(necking_region(:, 1), necking_region(:, 2), 'r');
plot(post_UTS_region(:, 1), post_UTS_region(:, 2), 'r'); % Plot post-UTS points as necking region (red)
xlabel('Engineering Strain');
ylabel('Engineering Stress');
legend('Compliance', 'Elastic Deformation', 'Necking');

% Add visual markers and tie-lines for UTS
plot(strain(idx_UTS), UTS, 'kx', 'MarkerSize', 10, 'LineWidth', 2);
plot([strain(idx_UTS), strain(idx_UTS)], [0, UTS], 'k--');

% Add labels for the beginning and end of each region
text(compliance_region(1, 1), compliance_region(1, 2), 'Compliance Start', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
text(compliance_region(end, 1), compliance_region(end, 2), 'Compliance End', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
text(elastic_region(1, 1), elastic_region(1, 2), 'Elastic Start', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
text(elastic_region(end, 1), elastic_region(end, 2), 'Elastic End', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
text(necking_region(1, 1), necking_region(1, 2), 'Necking Start', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
text(necking_region(end, 1), necking_region(end, 2), 'Necking End', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
text(strain(idx_UTS), UTS, sprintf(' UTS (%.2f, %.2f)', strain(idx_UTS), UTS), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

% Calculate the x-coordinate of the ideal line's starting point on the x-axis (y = 0)
x_ideal_start = elastic_region(1, 1) - (elastic_region(1, 2) / elastic_slope);

% Plot the ideal start of the test from the x-axis to the start of the elastic deformation region
x_ideal_line = linspace(x_ideal_start, elastic_region(1, 1), 100);
y_ideal_line = elastic_slope * (x_ideal_line - x_ideal_start);
plot(x_ideal_line, y_ideal_line, 'g--');



hold off;

% Set the size of the tiled layout window
figure('Position', [100, 100, 1800, 600]);

% Create a tiled layout for the three plots
tiledlayout(1, 3);

% Plot the original stress-strain curve in the first tile
nexttile;
hold on;
plot(compliance_region(:, 1), compliance_region(:, 2), 'b');
plot(elastic_region(:, 1), elastic_region(:, 2), 'g');
plot(necking_region(:, 1), necking_region(:, 2), 'r');
plot(post_UTS_region(:, 1), post_UTS_region(:, 2), 'r');
plot(x_ideal_line, y_ideal_line, 'g--');
xlabel('Engineering Strain');
ylabel('Engineering Stress');
title('Original Stress-Strain Curve');
xlim([0, max(all_files_data.data_5N_pure_10.Engr_Strain)]);
hold off;

% Create a new array of values for strain called offset_ideal_strain_
offset_ideal_strain_ = all_files_data.data_5N_pure_10.Engr_Strain - x_ideal_start;

% Offset the x-coordinates of the ideal line
x_ideal_line_offset = x_ideal_line - x_ideal_start;

% Create the new plot with offset strain values in the second tile
nexttile;
hold on;
plot(offset_ideal_strain_(1:size(compliance_region, 1)), compliance_region(:, 2), 'b');
plot(offset_ideal_strain_(size(compliance_region, 1) + 1:size(compliance_region, 1) + size(elastic_region, 1)), elastic_region(:, 2), 'g');
plot(offset_ideal_strain_(size(compliance_region, 1) + size(elastic_region, 1) + 1:size(compliance_region, 1) + size(elastic_region, 1) + size(necking_region, 1)), necking_region(:, 2), 'r');
plot(offset_ideal_strain_(size(compliance_region, 1) + size(elastic_region, 1) + size(necking_region, 1) + 1:end), post_UTS_region(:, 2), 'r');
plot(x_ideal_line_offset, y_ideal_line, 'g--');
xlabel('Offset Engineering Strain');
ylabel('Engineering Stress');
title('Stress-Strain Curve with Offset Strain');
xlim([0, max(offset_ideal_strain_)]);
legend('Compliance', 'Elastic Deformation', 'Necking');
hold off;


% Create the third plot with smoothed, offset stress-strain data
nexttile;
hold on;

% Smooth stress data and plot
window_size = 100;
smooth_compliance_stress = smoothdata(compliance_region(:, 2), 'movmean', window_size);
smooth_elastic_stress = smoothdata(elastic_region(:, 2), 'movmean', window_size);
smooth_necking_stress = smoothdata(necking_region(:, 2), 'movmean', window_size);
smooth_post_UTS_stress = smoothdata(post_UTS_region(:, 2), 'movmean', window_size);

% Concatenate smoothed stress data
smooth_stress = [smooth_compliance_stress; smooth_elastic_stress; smooth_necking_stress; smooth_post_UTS_stress];
smooth_strain = [compliance_region(:, 1); elastic_region(:, 1); necking_region(:, 1); post_UTS_region(:, 1)];

% Offset strain data
offset_smooth_strain = smooth_strain - x_ideal_start;

% Plot smoothed, offset stress-strain curve with coloring and increased line weight
plot(offset_smooth_strain(1:length(smooth_compliance_stress)), smooth_compliance_stress, 'b', 'LineWidth', 2);
plot(offset_smooth_strain(length(smooth_compliance_stress)+1:length(smooth_compliance_stress)+length(smooth_elastic_stress)), smooth_elastic_stress, 'g', 'LineWidth', 2);
plot(offset_smooth_strain(length(smooth_compliance_stress)+length(smooth_elastic_stress)+1:end), [smooth_necking_stress; smooth_post_UTS_stress], 'r', 'LineWidth', 2);
plot(x_ideal_line - x_ideal_start, y_ideal_line, 'g--', 'LineWidth', 2);

% Fill gaps between colored sections
line([offset_smooth_strain(length(smooth_compliance_stress)), offset_smooth_strain(length(smooth_compliance_stress)+1)], ...
    [smooth_compliance_stress(end), smooth_elastic_stress(1)], 'Color', 'b', 'LineWidth', 2);
line([offset_smooth_strain(length(smooth_compliance_stress)+length(smooth_elastic_stress)), offset_smooth_strain(length(smooth_compliance_stress)+length(smooth_elastic_stress)+1)], ...
    [smooth_elastic_stress(end), smooth_necking_stress(1)], 'Color', 'g', 'LineWidth', 2);


xlabel('Engineering Strain (Offset)');
ylabel('Engineering Stress');
title('Smoothed, Offset Stress-Strain Curve');
xlim([0, max(offset_smooth_strain)]);

xlabel('Engineering Strain (Offset)');
ylabel('Engineering Stress');
title('Smoothed, Offset Stress-Strain Curve');
xlim([0, max(offset_smooth_strain)]);


% Calculate R-squared values and standard deviations
Rsq_compliance = 1 - sum((compliance_region(:, 2) - smooth_compliance_stress).^2) / sum((compliance_region(:, 2) - mean(compliance_region(:, 2))).^2);
Rsq_elastic = 1 - sum((elastic_region(:, 2) - smooth_elastic_stress).^2) / sum((elastic_region(:, 2) - mean(elastic_region(:, 2))).^2);
Rsq_necking = 1 - sum((necking_region(:, 2) - smooth_necking_stress).^2) / sum((necking_region(:, 2) - mean(necking_region(:, 2))).^2);

std_compliance = std(compliance_region(:, 2));
std_elastic = std(elastic_region(:, 2));
std_necking = std(necking_region(:, 2));

% Add an information box with fit details
annotation('textbox', [0.7, 0.2, 0.2, 0.3], 'String', ...
    {sprintf('Compliance Region: R^2 = %.4f, Std = %.4f', Rsq_compliance, std_compliance), ...
     sprintf('Elastic Region: R^2 = %.4f, Std = %.4f', Rsq_elastic, std_elastic), ...
     sprintf('Necking Region: R^2 = %.4f, Std = %.4f', Rsq_necking, std_necking)});

hold off;




%**Next Section



%**************************************************************************************
% Function to toggle the visibility of a plot when the associated checkbox is changed
function toggle_plot(src, ~, plot_handle)
    if src.Value
        set(plot_handle, 'Visible', 'on');
    else
        set(plot_handle, 'Visible', 'off');
    end
end


