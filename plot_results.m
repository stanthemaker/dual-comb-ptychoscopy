function plot_results(file_list)
% plot_reconstruction_results Plots multiple reconstruction .mat files.
% The legend labels are automatically extracted from the file names.
%
% INPUTS:
%   file_list - Cell array of strings containing the file paths to load.
%
% EXAMPLE USAGE:
%   files = {'simulation_results/result_k=5.mat', 'simulation_results/result_k=20.mat'};
%   plot_reconstruction_results(files);

    num_files = length(file_list);

    % --- Setup Plot ---
    figure('Position', [100, 100, 800, 500]);
    hold on;
    
    % Get a color palette for an arbitrary number of lines
    colors = lines(num_files); 

    % --- Load Ground Truth ---
    try
        first_data = load(file_list{1});
        freq_THz = first_data.global_freq_axis / 1e12;
        
        % Plot ground truth (Black dashed line)
        plot(freq_THz, first_data.S_ss_global, 'r--', 'LineWidth', 2, ...
            'DisplayName', 'Ground Truth (S_{ss})');
    catch
        error('Failed to load the first file. Please check the file path: %s', file_list{1});
    end

    % --- Loop through and plot all files ---
    for i = 1:num_files
        try
            data = load(file_list{i});
            
            % Extract just the file name (and extension) from the full path
            [~, fileName, ext] = fileparts(file_list{i});
            label_name = [fileName, ext]; % e.g., 'result_k=5.mat'
            % Note: If you want to drop the '.mat' in the legend, just use: label_name = fileName;
            
            plot(freq_THz, data.spec_recon_dB, 'Color', colors(i,:), ...
                'LineWidth', 1, 'DisplayName', label_name);
        catch
            warning('Could not load file: %s. Skipping...', file_list{i});
        end
    end

    % --- Formatting --
    xlabel('Optical Frequency (THz)');
    ylabel('Signal Power (dB)');
    xlim([192.7 193.4]);
    ylim([-30 10]); 
    grid on;
    
    % Use Interpreter 'none' so underscores in filenames don't become subscripts!
    legend('Location', 'best', 'Interpreter', 'none');
    set(gca, 'FontSize', 14);
    
    hold off;
end