close all;
clear;


load Plot_Data_For_APS_Paper.mat % load data 
addpath("utils_Ryan\misc_utils\");
%% Plot Processed B-Scans

Final_Pulsed_x = normalize_to_range(Final_Pulsed_x, 0 , 1.4);
Final_FMCW_x = normalize_to_range(Final_FMCW_x, 0 , 1.4);


figure();
plot_D(normalize(Final_FMCW_Processed_B_Scan), Final_FMCW_x,Final_FMCW_t.*10^18,"")
set(gca,'YDir','reverse');
xlabel("$$\bf{Scan\: x \:(m)}$$" , 'FontWeight','bold' ,FontSize=20, FontName="times");
ylabel("$$\bf{}$$", 'FontWeight','bold',FontSize=20, FontName="times");
set(gca,'FontWeight','bold', FontName="times")
xticks([0 0.7 1.4])
yticks([])
colorbar;
%xticks([0 0.5 1 1.3])

final_CFAR = cfar_detector(Final_FMCW_Processed_B_Scan, 20, 2, 30000);
final_CFAR = final_CFAR - mean(final_CFAR,2);

figure();
plot_D(normalize(final_CFAR), Final_FMCW_x,Final_FMCW_t.*10^18,"")
set(gca,'YDir','reverse');
xlabel("$$\bf{Scan\: x \:(m)}$$" , 'FontWeight','bold' ,FontSize=20, FontName="times");
ylabel("$$\bf{}$$", 'FontWeight','bold',FontSize=20, FontName="times");
set(gca,'FontWeight','bold', FontName="times")
xticks([0 0.7 1.4])
yticks([])

%% Plotting Functions Change these

function Mixed_cfar = cfar_detector(Mixed, N_train, N_guard, P_fa)
    % CFAR Detector for Frequency Spectrum (Column-wise Processing)
    % Inputs:
    %   - Mixed: Input matrix (each column processed separately)
    %   - N_train: Number of training cells
    %   - N_guard: Number of guard cells
    %   - P_fa: False alarm probability
    % Output:
    %   - Mixed_cfar: CFAR-detected matrix (same size as Mixed)
    % Normalize input signal (column-wise)
    Mixed = abs(Mixed) ./ max(abs(Mixed), [], 1); % Normalize each column
    % Compute CFAR threshold factor
    alpha = (N_train * (P_fa^(-1/N_train)) - 1);
    % Initialize CFAR output matrix
    Mixed_cfar = zeros(size(Mixed));
    % Get matrix dimensions
    [num_rows, num_cols] = size(Mixed);
    % Process each column independently
    for col = 1:num_cols
        for i = (N_train + N_guard + 1):(num_rows - (N_train + N_guard))
            % Training region (excluding guard cells)
            training_cells = [Mixed(i - N_train - N_guard : i - N_guard - 1, col); ...
                              Mixed(i + N_guard + 1 : i + N_train + N_guard, col)];
            % Estimate noise level using mean
            noise_power = mean(training_cells);
            % Compute CFAR threshold
            threshold = alpha * noise_power;
            % Apply CFAR detection
            if Mixed(i, col) > threshold
                Mixed_cfar(i, col) = Mixed(i, col); % Keep detected target
            end
        end
    end
end








function normalized_vector = normalize_to_range(input_vector, new_min, new_max)
    % Normalize a vector to a specified range [new_min, new_max]
    % 
    % Parameters:
    %   input_vector - The vector to be normalized
    %   new_min - The minimum value of the new range
    %   new_max - The maximum value of the new range
    %
    % Returns:
    %   normalized_vector - The input vector normalized to the specified range
    
    % Find the current range of the input vector
    old_min = min(input_vector);
    old_max = max(input_vector);
    
    % Scale and shift the input vector
    normalized_vector = ((input_vector - old_min) / (old_max - old_min)) * ...
                         (new_max - new_min) + new_min;
end


function x = normalize(x)
    x = x./max(abs(x(:)));
end

function [ax] = plot_I(I,x,z,DR,title_str)
colormap(jet);

imagesc(x,z.*1e2,10.*log10(abs(normalize(I))));
ax = gca();
xlabel("x (m)");
ylabel("z (cm)");
title(title_str);
caxis([-DR 0]);


end

function [ax] = plot_D(D,x,t,title_str)
    % create colorbar limits
    clim_abs_max = 250;
    threshold    = 0.99;
    D_sort = sort(D(:),'ascend');
    clims  = [D_sort(round(numel(D_sort)*(1-threshold))), D_sort(round(numel(D_sort)*threshold))];
    clims  = [max(-clim_abs_max, clims(1)),min(clim_abs_max, clims(2))];

    colormap(gray(2^12));
    imagesc(x, t/1e9, D);
    ax = gca();
    %colorbar;
    caxis(clims);
    title(title_str, FontSize=20);
    xlabel("x (m)", FontSize=20);
    ylabel("f (Ghz)", FontSize=20);
end
